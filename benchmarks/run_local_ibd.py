#!/usr/bin/env python3
"""Run autoreporting locally on a subset (one chromosome / region) of an external sumstat,
mirroring the WDL command in autoreporting_external.json, with timing + cProfile.

Why: the WDL `report` task is slow on the many-variant/many-signal external sumstats and
the tabix LD step dominates. This driver reproduces a single phenotype run locally so the
LD/grouping cost can be profiled and iterated on, without spinning up Cromwell or the docker.

It resolves inputs from the WDL json (+ the input array tsv), tabix-slices the requested
chromosome straight from GCS (local pysam here has htslib GCS support), runs Scripts/main.py
on the slice (LD panel stays remote), and prints a pstats summary when --profile is set.

Examples:
  # profile the LD/grouping path on chr21, serial (cProfile sees only the main process)
  python3.10 benchmarks/run_local_ibd.py --chrom 21 --profile

  # measure parallel wall-clock + the variant1-indexed narrow-fetch win
  python3.10 benchmarks/run_local_ibd.py --chrom 21 --workers 8 --assume-variant1-indexed

  # full pipeline (annotations + gwas catalog), small region
  python3.10 benchmarks/run_local_ibd.py --chrom 21 --region 5000000-15000000 --mode full
"""
import argparse
import json
import os
import subprocess
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
MAIN = os.path.join(REPO, "Scripts", "main.py")
DEFAULT_JSON = "/home/jkarjala/suite/genetics-results-munge/wdl/autoreporting_external.json"


def gcs_token():
    return subprocess.run(["gcloud", "auth", "print-access-token"],
                          capture_output=True, encoding="utf-8", check=True).stdout.strip()


def resolve_local(gs_path, json_dir):
    """Prefer a local copy sitting next to the json (the repo keeps the input tsv / pheno
    info there); otherwise return the gs:// path for `gsutil cat`."""
    local = os.path.join(json_dir, os.path.basename(gs_path))
    return local if os.path.exists(local) else gs_path


def read_input_row(input_tsv, pheno, json_dir):
    path = resolve_local(input_tsv, json_dir)
    if path.startswith("gs://"):
        text = subprocess.run(["gsutil", "cat", path], capture_output=True,
                              encoding="utf-8", check=True).stdout
        lines = text.splitlines()
    else:
        with open(path) as f:
            lines = f.read().splitlines()
    for line in lines:
        parts = line.split("\t")
        if parts and parts[0] == pheno:
            return parts
    raise SystemExit(f"phenotype {pheno!r} not found in {path}")


def build_subset(sumstat_gs, contig, region, outdir, rebuild):
    """Tabix-slice `contig` (optionally a start-end region) from the remote sumstat into a
    local bgzipped+indexed tsv. Cached by name; pass --rebuild-subset to force."""
    tag = contig if region is None else f"{contig}_{region[0]}_{region[1]}"
    base = os.path.join(outdir, f"subset_{tag}.tsv")
    gz = base + ".gz"
    if os.path.exists(gz) and os.path.exists(gz + ".tbi") and not rebuild:
        print(f"[subset] reusing cached {gz}")
        return gz
    print(f"[subset] slicing {sumstat_gs} contig={contig} region={region} ...", flush=True)
    tf = pysam.TabixFile(sumstat_gs, encoding="utf-8")
    header = tf.header[-1] if tf.header else None
    if header is None:
        raise SystemExit("source sumstat has no embedded header")
    start, end = (None, None) if region is None else region
    n = 0
    t0 = time.perf_counter()
    with open(base, "w") as out:
        out.write(header + "\n")
        for row in tf.fetch(contig, start, end):
            out.write(str(row) + "\n")
            n += 1
    print(f"[subset] wrote {n} rows in {time.perf_counter()-t0:.1f}s, indexing ...", flush=True)
    # sumstat layout: chrom=col0, pos=col1; header line starts with '#'
    pysam.tabix_index(base, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    return gz


def build_argv(cfg, row, subset_gz, out_prefix, args):
    rep = lambda k: cfg[f"autoreporting.report.{k}"]
    ld_panel = cfg["autoreporting.ld_panel"]
    ld_api = cfg["autoreporting.ld_api"]
    locus_width_kb = cfg["autoreporting.locus_width_kb"]
    column_names = rep("column_names")
    pheno = row[0]

    argv = [
        subset_gz,
        "--pheno-name", pheno,
        "--sign-treshold", str(rep("sign_treshold")),
        "--alt-sign-treshold", str(rep("alt_sign_treshold")),
        "--group",
        "--grouping-method", "ld",
        "--locus-width-kb", str(locus_width_kb),
        "--ld-panel-path", ld_panel,
        "--ld-api", ld_api,
        "--plink-memory", str(rep("plink_memory")),
        "--strict-group-r2", str(rep("strict_group_r2")),
        "--column-labels", *column_names,
        "--extra-cols", *rep("extra_columns").split(),
        "--report-out", out_prefix + ".report.out",
        "--top-report-out", out_prefix + ".top.out",
        "--ld-workers", str(args.workers),
    ]
    # ld_opts from the json, e.g. "--dynamic-r2-chisq 5.0"
    argv += rep("ld_opts").split()
    if rep("pval_is_mlog10p"):
        argv.append("--pval-is-mlog10p")
    if rep("overlap"):
        argv.append("--overlap")
    # narrow fetch is on by default in main.py now; pass the explicit opt-out to benchmark
    # the old wide path
    argv.append("--ld-assume-variant1-indexed" if args.assume_variant1_indexed
                else "--no-ld-assume-variant1-indexed")

    json_dir = os.path.dirname(os.path.abspath(args.json))
    pheno_info = resolve_local(rep("phenotype_info_file"), json_dir)
    if not pheno_info.startswith("gs://"):
        argv += ["--pheno-info-file", pheno_info]

    if args.mode == "full":
        argv += [
            "--finngen-path", rep("finngen_annotation"),
            "--functional-path", rep("functional_annotation"),
            "--gnomad-path", cfg["autoreporting.report.gnomad"],
            "--use-gwascatalog",
            "--gwascatalog-threads", str(rep("gwascatalog_threads")),
            "--gwascatalog-pval", str(rep("gwascatalog_pval")),
            "--gwascatalog-width-kb", str(rep("gwascatalog_width_kb")),
            "--db", rep("db_choice"),
            "--local-gwascatalog", rep("local_gwcatalog"),
            "--gwascatalog-allele-file", rep("allele_vcf_file"),
            "--custom-dataresource", rep("custom_dataresource"),
        ]
        if rep("finngen_variants_only"):
            argv.append("--finngen-variants-only")
    return argv


def print_profile(prof_path, top):
    import pstats
    print("\n" + "=" * 80 + f"\nPROFILE (top {top} by cumulative time)\n" + "=" * 80)
    st = pstats.Stats(prof_path)
    st.strip_dirs().sort_stats("cumulative").print_stats(top)
    print("=" * 80 + f"\nPROFILE (top {top} by total/own time)\n" + "=" * 80)
    st.sort_stats("tottime").print_stats(top)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", default=DEFAULT_JSON, help="WDL input json")
    ap.add_argument("--pheno", default="IBD", help="phenotype name (row in input array tsv)")
    ap.add_argument("--chrom", default="21", help="contig to slice (as named in the sumstat)")
    ap.add_argument("--region", default=None, help="optional START-END within the chrom")
    ap.add_argument("--mode", choices=["group", "full"], default="group",
                    help="group: parse+LD grouping+report only (LD focus); full: + annotations & gwas catalog")
    ap.add_argument("--workers", type=int, default=1, help="--ld-workers (1 = serial; use 1 with --profile)")
    ap.add_argument("--assume-variant1-indexed", action=argparse.BooleanOptionalAction, default=True,
                    help="1bp narrow LD fetch (valid for finngen panels). On by default; "
                         "--no-assume-variant1-indexed benchmarks the old wide ±range fetch")
    ap.add_argument("--profile", action="store_true", help="run under cProfile and print a pstats summary")
    ap.add_argument("--top", type=int, default=30, help="rows in the pstats summary")
    ap.add_argument("--outdir", default=os.path.join(HERE, "_out"))
    ap.add_argument("--rebuild-subset", action="store_true")
    args = ap.parse_args()

    if args.profile and args.workers != 1:
        print("WARNING: cProfile only sees the main process; LD worker time will hide behind "
              "the prefetch wait. Use --workers 1 for an accurate function breakdown.",
              file=sys.stderr)

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    # work from outdir so pysam's remote .tbi index caches land here (gitignored) instead of
    # polluting the repo root. all other paths used below are absolute.
    os.chdir(outdir)
    args.outdir = outdir
    with open(args.json) as f:
        cfg = json.load(f)
    json_dir = os.path.dirname(os.path.abspath(args.json))
    row = read_input_row(cfg["autoreporting.input_array_file"], args.pheno, json_dir)
    sumstat_gs = row[1]

    region = None
    if args.region:
        lo, hi = args.region.split("-")
        region = (int(lo), int(hi))

    # set on this process's env too: build_subset opens gs:// via pysam here, and the token
    # is inherited by the main.py subprocess
    os.environ["GCS_OAUTH_TOKEN"] = gcs_token()
    env = dict(os.environ)

    subset_gz = build_subset(sumstat_gs, args.chrom, region, args.outdir, args.rebuild_subset)

    out_prefix = os.path.join(args.outdir, f"{args.pheno}_{args.mode}")
    argv = build_argv(cfg, row, subset_gz, out_prefix, args)

    prof_path = out_prefix + ".prof"
    cmd = [sys.executable]
    if args.profile:
        cmd += ["-m", "cProfile", "-o", prof_path]
    cmd += [MAIN] + argv

    print("\n[run] mode={} workers={} assume_variant1_indexed={} profile={}".format(
        args.mode, args.workers, args.assume_variant1_indexed, args.profile))
    print("[run] " + " ".join(cmd) + "\n", flush=True)

    t0 = time.perf_counter()
    # cwd is outdir (set above) so remote .tbi caches land there; main.py is invoked by
    # absolute path, so Scripts/ is still added to sys.path[0] for its imports
    rc = subprocess.run(cmd, env=env).returncode
    wall = time.perf_counter() - t0
    print(f"\n[run] return code {rc}, total wall time {wall:.1f}s")

    if args.profile and rc == 0:
        print_profile(prof_path, args.top)
    if rc != 0:
        sys.exit(rc)


if __name__ == "__main__":
    main()
