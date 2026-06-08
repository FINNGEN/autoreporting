# Local benchmarking / profiling

`run_local_ibd.py` reproduces a single autoreporting `report`-task run locally on a subset
of an external sumstat (one chromosome or region), mirroring the command in
`wdl/autoreporting_external.json`. Use it to profile and iterate on the LD/grouping cost
without Cromwell or docker.

It resolves inputs from the WDL json (+ the input array tsv), tabix-slices the chosen
chromosome straight from GCS (local pysam/htslib here has GCS support), runs
`Scripts/main.py` on the slice (the LD panel stays remote), and prints a `pstats` summary
when `--profile` is set. Stage timings are logged with a greppable `[STAGE]` prefix from
`time_decorator.timed`.

## Usage

```bash
# profile the LD/grouping path on chr21 (serial: cProfile only sees the main process)
python3.10 benchmarks/run_local_ibd.py --chrom 21 --profile

# small region, narrow LD fetch + parallel workers (wall-clock)
python3.10 benchmarks/run_local_ibd.py --chrom 21 --region 5000000-15000000 \
    --workers 8 --assume-variant1-indexed

# full pipeline (annotations + gwas catalog)
python3.10 benchmarks/run_local_ibd.py --chrom 21 --region 5000000-15000000 --mode full
```

Flags: `--mode group|full` (default `group` = parse + LD grouping + report only),
`--workers N` (`--ld-workers`; use `1` with `--profile`), `--region START-END`, `--pheno`,
`--rebuild-subset`. The narrow 1bp LD fetch is on by default (matching the code default);
`--no-assume-variant1-indexed` benchmarks the old wide ±range fetch. Outputs and the cached
slice land in `benchmarks/_out/` (gitignored).

## Finding (IBD, chr21:5,000,000-15,000,000, 5 GWS leads, serial)

| LD fetch mode | LD prefetch time | total wall |
|---|---|---|
| wide ±2Mb window (default)            | **88.5 s** | 95.5 s |
| narrow 1bp (`--assume-variant1-indexed`) | **0.29 s** | 3.9 s |

Reports were **byte-identical** between the two modes on the real finngen_r12 panel.

The cProfile of the wide path shows why: fetching each lead's ±2Mb window from the 1.6 GB
LD file iterates **~94 million tabix rows across just 5 leads** —
`libctabix __next__/__cnext__` accounts for ~104 s and the Python parse in
`TabixLD.get_range` ~26 s. The narrow fetch eliminates the scan because the panel is indexed
by `variant1`'s position (verified: the file's `pos` column equals variant1's position in
every row), so all of a lead's rows sit at a single coordinate.

The narrow fetch is now the **default** (`assume_variant1_indexed=True`) since it is the
dominant win and verified correct on the real finngen_r12 panel. Pass
`--no-ld-assume-variant1-indexed` (main.py / post_process_hits.py) only for a panel not
indexed by variant1 position. The other optimizations (asTuple parsing, header-index
hoisting, greedy de-quadratic) help the parse/grouping stages and the wide path.
