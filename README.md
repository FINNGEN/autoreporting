# Automatic hit reporting tool for FINNGEN

## A tool for finding gws results from GWAS summary statistics

This pipeline is used to
1) Filter out and group gws variants from FINNGEN summary statistics 
2) Perform finemapping on the filtered SNPs (TBD)
3) Annotate the gws variants using gnoMAD and FINNGEN annotations 
4) Cross-reference the variants to previous results, e.g. gwascatalog summary statistic database or hand-picked results from studies 
Currently, steps 1,3 and 4 are operational. 

## Dependencies

packages:python 3, plink 1.9, ldstore (tested on 1.1)

python 3 libraries: requests, numpy, pandas, pytabix 

## Installation

Install dependencies

```
sudo apt install tabix python3 pip plink1.9
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz
tar xvf ldstore_v1.1_x86_64.tgz
sudo cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/ldstore
pip3 install pytabix requests numpy pandas
```

Copy repository to folder

```
git clone https://github.com/FINNGEN/autoreporting.git
```

## Usage

In the project folder, the script can be used by either calling the whole script or calling the individual scripts by themselves.
### main:

```
usage: main.py [-h] [--sign-treshold SIG_TRESHOLD] [--fetch_out FETCH_OUT]
               [--group] [--grouping-method GROUPING_METHOD]
               [--locus-width-kb LOC_WIDTH]
               [--alt-sign-treshold SIG_TRESHOLD_2]
               [--ld-panel-path LD_PANEL_PATH] [--ld-r2 LD_R2]
               [--plink-memory PLINK_MEM] [--overlap]
               [--gnomad-genome-path GNOMAD_GENOME_PATH]
               [--gnomad-exome-path GNOMAD_EXOME_PATH] [--include-batch-freq]
               [--finngen-path FINNGEN_PATH] [--annotate-out ANNOTATE_OUT]
               [--compare-style COMPARE_STYLE]
               [--summary-fpath FILE [FILE ...]]
               [--endpoints ENDPOINTS [ENDPOINTS ...]] [--build-38]
               [--check-for-ld] [--raport-out RAPORT_OUT]
               [--gwascatalog-pval GWASCATALOG_PVAL]
               [--gwascatalog-width-kb GWASCATALOG_PAD]
               [--ldstore-threads LDSTORE_THREADS] [--ld-treshold LD_TRESHOLD]
               gws_fpath
```

Command-line arguments:
Argument   |  Meaning   |   Example | Original script
--- | --- | --- | ---
--sign-treshold | signifigance treshold for variants | --sign-treshold 5e-8 | gws_fetch.py
--fetch-out | output file path for filtered and/or grouped variants. 'fetch_out.csv' by default. | --fetch-out output.tsv | gws_fetch.py
--group | supplying this flag results in the variants being grouped into groups | --group | gws_fetch.py
--grouping-method | grouping method used if --group flag is supplied. options are 'simple', i.e. grouping based on range from gws variants, or 'ld', i.e. grouping using plink --clump |  --grouping-method ld | gws_fetch.py
--locus-width-kb | group widths in kb. In case of ld clumping, the value is supplied to plink --clump-kb. | --locus-width-kb 500 | gws_fetch.py
--alt-sign-treshold | optional alternate signifigance treshold for including less significant variants into groups. | --alt-sign-treshold 5e-6 | gws_fetch.py
--ld-panel-path | path to ld panel, without panel file suffix. Ld panel must currently be in plink's .bed format, as a single file. | --ld-panel-path path_to_panel/plink_file | gws_fetch.py
--ld-r2 | plink clump-r2 argument, default 0.4 | --ld-r2 0.7 | gws_fetch.py
--plink-memory | plink --memory argument. Default 12000 | --plink-memory 16000 | gws_fetch.py
--overlap | If this flag is supplied, the groups of gws variants are allowed to overlap, i.e. a single variant can appear multiple times in different groups. | --overlap | gws_fetch.py
--gnomad-genome-path | path to gnomad genome annotation file. Must be tabixed. Required for annotation. | --gnomad-genome-path gnomad_path/gnomad_file.tsv.gz | annotate<span></span>.py
--gnomad-exome-path | path to gnomad exome annotation file. Must be tabixed. Required for annotation. | --gnomad-exome-path gnomad_path/gnomad_file.tsv.gz | annotate<span></span>.py
--include-batch-freq | Include batch frequencies from finngen annotation file | --include-batch-freq | annotate<span></span>.py
--finngen-path | Path to finngen annotation file, containing e.g. most severe consequence and corresponding gene of the variants | --finngen-path path_to_file/annotation.tsv.gz | annotate<span></span>.py
--annotate-out | annotation output file, default 'annotate_out.csv' | --annotate-out annotation_output.tsv | annotate<span></span>.py
--compare-style | Whether to use gwascatalog or additional summary statistics to compare findings to literature. Use values 'file' or 'gwascatalog', default 'gwascatalog' | --compare-style 'gwascatalog' | compare<span></span>.py
--summary-fpath | filepaths to external summary statistic files. Must be in GRCh38. Must be same amount as endpoints. Do not supply gws_path directly after this one. | --summary-fpath first_summary.tsv second_summary.tsv third_summary.tsv | compare<span></span>.py
--endpoints | Phenotypes/endpoints for the supplied summary statistics. Must be in same order as summary statistics. Do not supply gws_path directly after this one. | --endpoints first_endpoint second_endpoint third_endpoint | compare<span></span>.py
--build-38 | Supply if summary statistics are in GRCh38. Currently required for external summary statistics. | --build-38 | compare<span></span>.py
--check-for-ld | When supplied, gws variants and summary statistics (from file or gwascatalog) are tested for ld using LDstore.  | --check-for-ld | compare<span></span>.py
--gwascatalog-pval | P-value to use for filtering results from gwascatalog's summary statistic API. default 5e-8 | --gwascatalog-pval 5e-6 | compare<span></span>.py
--gwascatalog-width-kb | Buffer outside gws variants that is searched from gwascatalog, in kb. Default 25  | --gwascatalog-width-kb 50 | compare<span></span>.py
--ldstore-threads | Number of threads to use with LDstore. Currently not used. | --ldstore-threads 1 | compare<span></span>.py
--ld-treshold | LD treshold for LDstore, above of which summary statistic variants in ld with our variants are included. Default 0.4 | --ld-treshold 0.8 | compare<span></span>.py
gws_path |  Path to the tabixed and gzipped summary statistic that is going to be filtered, annotated and compared. Required argument. | path_to_summary_statistic/summary_statistic.tsv.gz | gws_fetch.py

The same arguments are used in the smaller scripts that the main script uses.

For example, In order to filter genome-wide significant variants, and to compare them against gwascatalog's summary statistics API, the following call can be used:
```
python3 Scripts/main.py path_to_ss/summ_stat.tsv.gz
```

Additional features, such as result grouping, can be added through the use of the aforementioned flags.

### gws_fetch.py:

```
usage: gws_fetch.py [-h] [--sign-treshold SIG_TRESHOLD]
                    [--fetch_out FETCH_OUT] [--group]
                    [--grouping-method GROUPING_METHOD]
                    [--locus-width-kb LOC_WIDTH]
                    [--alt-sign-treshold SIG_TRESHOLD_2]
                    [--ld-panel-path LD_PANEL_PATH] [--ld-r2 LD_R2]
                    [--plink-memory PLINK_MEM] [--overlap]
                    gws_fpath
```
 The gws_fetch.py-script is used to filter genome-wide significant variants from the summary statistic file as well as optionally group the variants, either based on a range around top hits, or by using plink's --clump functionality. The arguments used are the same as the ones in main<span></span>.py. For example, to just filter variants according to a p-value, the script can be called by
 ```
 python3 gws_fetch.py --sign-treshold 2.5e-8 path_to_ss/summary_statistic.tsv.gz
 ```

 ### annotate<span></span>.py:

```
usage: annotate.py [-h] [--gnomad-genome-path GNOMAD_GENOME_PATH]
                   [--gnomad-exome-path GNOMAD_EXOME_PATH]
                   [--include-batch-freq] [--finngen-path FINNGEN_PATH]
                   [--annotate-out ANNOTATE_OUT]
                   annotate_fpath
```

The annotate<span></span>.py-script is used to annotate the previously filtered genome-wide significant variants, using annotation files from gnoMAD as well as annotation files specific to the FINNGEN project. The arguments are the same as in main<span></span>.py, except for annotate_fpath, which is the path to the filtered variants. For example, to annotate variants the script can be called like this:
```
python3 annotate.py variant_file_path/variants.tsv --gnomad-genome-path path_to_gnomad/gnomad_genomes.tsv.gz --gnomad-exome-path path_to_gnomad/gnomad_exomes.tsv.gz --finngen-path path_to_finngen_annotation/finngen_annotation.tsv.gz --annotate-out annotation_output.tsv
```

### compare<span></span>.py:

```
usage: compare.py [-h] [--compare-style COMPARE_STYLE]
                  [--summary-fpath FILE [FILE ...]]
                  [--endpoints ENDPOINTS [ENDPOINTS ...]] [--build-38]
                  [--check-for-ld] [--ld-panel-path LD_PANEL_PATH]
                  [--raport-out RAPORT_OUT]
                  [--gwascatalog-pval GWASCATALOG_PVAL]
                  [--gwascatalog-width-kb GWASCATALOG_PAD]
                  [--ldstore-threads LDSTORE_THREADS]
                  [--ld-treshold LD_TRESHOLD]
                  compare_fname
```
The compare<span></span>.py-script is used to compare the genome-wide significant variants to earlier results, either in the form of summary statistics supplied to the script or searched from GWAScatalog's summary statistic api. The arguments are the same as in main<span></span>.py, except for compare_fname, which is the input variant file. For example, to simply check if the variants have any corresponding hits in GWAScatalog summary statistics, one can use the following command:
```
python3 Scripts/compare.py variant_file.tsv --compare-style gwascatalog --gwascatalog-pval 5e-8 --raport-out output_raport.tsv
```
Additional flags, such as `--check-for-ld`, can be used to check if the summary statistics are in ld with the variants, and to report them if they are.



