# Automatic hit reporting tool for FINNGEN

## A tool for finding gws results from GWAS summary statistics

This pipeline is used to
1) Filter out and group gws SNPs from FINNGEN summary statistics
2) Perform finemapping on the filtered SNPs
3) Annotate the SNPs using gnoMAD and FINNGEN annotations
4) Cross-reference the SNPs to previous results, e.g. online databases and/or hand-picked results from studies

## 1: Filtering and clumping GWS variants: gws_fetch.py
The first part of the tool filters in genome-wide-significant variants (treshold configurable) and optionally groups them according to the command line settings. The grouping can be done using either grouping based on an area around top gws variants, or using PLINK's ld clumping. An alternate signifigance treshold for including variants in the groups can be specified.

### Usage

usage: gws_fetch.py [-h] [-s SIG_TRESHOLD] [-o OUT_FNAME] [-g]
                    [--grouping-method GROUPING_METHOD] [-w LOC_WIDTH]
                    [-s2 SIG_TRESHOLD_2] [--ld-panel-path LD_PANEL_PATH]
                    [--ld-r2 LD_R2] [--plink-memory PLINK_MEM]
                    gws_fpath
Argument list:
    gws_fpath   Filepath of the compressed tsv
    -h, --help            show this help message and exit
    --signifigance-treshold SIG_TRESHOLD    Signifigance treshold
    --out-fname OUT_FNAME   Output filename, default is out.csv
    --group Whether to group SNPs
    --grouping-method GROUPING_METHOD   Decide grouping method, simple or ld, default simple
    --locus-width-kb LOC_WIDTH  locus width to include for each SNP, in kb
    --alternate-sign-treshold SIG_TRESHOLD_2    optional group treshold
    --ld-panel-path LD_PANEL_PATH   Filename to the genotype data for ld calculation without suffix
    --ld-r2 LD_R2   r2 cutoff for ld clumping
    --plink-memory PLINK_MEM    plink memory for ld clumping, in MB

## 2: Finemapping the gws variants
TBD

## 3: Annotating the gws variants using gnoMAD annotations and FINNGEN annotations
The third part 

usage: annotate.py [-h] [--gnomad-genome-path GNOMAD_GENOME_PATH]
                   [--gnomad-exome-path GNOMAD_EXOME_PATH]
                   [--include-batch-freq] [--finngen-path FINNGEN_PATH]
                   [-o OUT_FNAME]
                   annotate_fpath

Annotate results using gnoMAD and additional annotations

Argument list:
    annotate_fpath        Filepath of the results to be annotated
    -h, --help            show this help message and exit
    --gnomad-genome-path GNOMAD_GENOME_PATH Gnomad genome annotation file filepath
    --gnomad-exome-path GNOMAD_EXOME_PATH   Gnomad exome annotation file filepath
    --include-batch-freq    Include batch frequencies from finngen annotations
    --finngen-path FINNGEN_PATH Finngen annotation file filepath
    --out-fname OUT_FNAME   Output filename, default is out.csv
