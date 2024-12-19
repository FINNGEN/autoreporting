# Automatic hit reporting tool for FinnGen


## Description

This pipeline is used to
1) Filter and group genome-wide significant variants from FinnGen summary statistics.
2) Annotate the genome-wide significant variants using gnomAD and FinnGen annotations, as well as GWAS Catalog and/or hand-curated annotations.

__NOTE: currently, only files which are in build 38 are supported. This concerns all of the input files.__

__NOTE: In case you do not have access to the Finnish LD panel (requires access to refinery), use the command-line option `--ld-api online` to circumvent that and still get access to Finnish LD.__ 

## Installation

### Dependencies

Packages:
```
python 3 (version 3.6+)
pip
plink (latest binary)
ldstore (tested on 1.1) (only on docker image, will be removed at some point)
zlib
```
Python 3 libraries:
``` 
requests
numpy
pandas 
pysam 
```

##  Installationprocess

Install dependencies (on ubuntu 18.04+):

```
git clone https://github.com/FinnGen/autoreporting.git
sudo apt install tabix python3 python3-pip zlib1g-dev
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip 
unzip plink_linux_x86_64_latest.zip 
sudo cp plink /usr/local/bin/
pip3 install -r autoreporting/requirements.txt
```

## Resources

Refer to the latest WDL json template for up-to-date resources.

__NOTE: DEPRECATED__
The resources to use this tool (gnomAD & FinnGen annotations, LD panel) can be found here:
- Summary statistics: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/summary_stats/release```
- Credible sets: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/finemap```
- gnomAD genome & exome annotations: ```gs://finngen-production-library-green/autoreporting_annotations/gnomad_data/```
- FinnGen annotations: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/annotations/```
- Functional annotations: ```gs://r4_data_west1/gnomad_functional_variants/fin_enriched_genomes_select_columns.txt.gz```
- LD panel (based on 1000 genomes data): ```gs://finngen-production-library-green/autoreporting_annotations/1kg_ld/whole_1k*```
- VCF allele panel for gwas catalog comparison: ```gs://finngen-production-library-green/autoreporting_annotations/vcf_allele_panel/dbsnp_hg38p13_b153_numerical.vcf.gz```


## Running the WDL for multiple endpoints

### Inputs
Required inputs:
- Summary statistics
Optional inputs:
- Previous release summary statistics (or any other summstat you want to compare to)
- Finemapping results from susie (filtered files preferable)
- Annotation files, if not up-to-date in latest WDL template

### Create WDL inputs
First, create lists of your summary statistics, and possible credible set .cred files, and possible previous release files.
If you do not have finemapping files or previous release files, create empty lists.
For example:
```sh
gsutil ls gs://summary-stat-path/files/*.gz > summ_stat_list
gsutil ls gs://credible-set-path/files/*.snp.filter.tsv > credible_set_files
gsutil ls gs://prev-release-summstat/files/*.gz > previous_release_files
```
If no credset and/or previous release files:
```sh
touch credible_set_files
touch previous_release_files
```

Also, prepare your dummy value for an empty file that will be localized in case the input is not available.
By default, you can use `gs://misc-analysis/EMPTY_FILE`.
### Create input array
Use Scripts/wdl_processing_scripts/pheno-credset_array.py:
```sh
python3 Scripts/wdl_processing_scripts/pheno-credset_array.py --phenotype-list summ_stat_list \
  --phenotype-prefix POSSIBLE_FILENAME_PREFIX \
  --credset-list credible_set_files \
  --credset-prefix "" \
  --prev-release previous_release_files \
  --prev-release-prefix POSSIBLE_PREVREL_PREFIX \
  --empty-file-path "gs://misc-analysis/EMPTY_FILE" \
  --out autoreporting_input_array.tsv \
  --only-cred # this filters the inputs to only those that have finemapping results. Useful in core analysis, but if there are no finemap results, leave it out 
```
In case you have finemapped the results and only want to analyse the results which had finemapping results, use the last argument `--only-cred`. If you did not finemap your results, do not include that flag.

### Fill in WDL parameters
Fill in correct parameters in autoreporting.template.json
- input array
- correct docker image, latest is in template
- annotation files, e.g. gnomad and variant annotations. Latest should be in template, but you can confirm from the maintainer.
- GWAS Catalog annotation files, like the gwas catalog release tsv and an allele vcf file
    - You can download the most recent GWAS Catalog annotation from https://www.ebi.ac.uk/gwas/docs/file-downloads (All associations v1.0.2)
- custom data resource, i.e. betamatching external annotation from https://github.com/FINNGEN/betamatch-ext-data
    - You can use e.g. https://github.com/FINNGEN/autoreporting/blob/master/custom_catalog_script/bmatch_to_singlefile.py to combine the data files to create the custom data resource file.
- correct LD panel, as of R10 sisuv4.2
- run parameters:
    - primary grouping method, most likely cred
    - significance thresholds, sign_threshold for group lead pval threshold (not in cred) and alt_sign_threshold for LD partner sign threshold
    - LD threshold and LD options, like whether to use a static LD threshold or dynamic r2 based on inv chisq of peak pval
    - grouping locus width, i.e. the window size (variants as far as that window in kb can be chosen, but not variants farther away than that)
    - strict group r2, for LD grouping top reports
    - whether to ignore a region or not. If finemapping already ignored a region (e.g. HLA), you might not have to. If using other grouping options, you might want to ignore HLA.
    - dummy file, used for missing data. Use the same value as you used in the array creation.


## Running the tool locally

You can run the autoreporting tool locally by installing the dependencies noted above.
Running the tool is as simple as 
```sh
python3 Scripts/main.py SUMMARY_STATISTIC [OPTIONS]
```
You can see all of the available options with `--help`:
```sh
python3 Scripts/main.py --help
```

###  Command-line arguments
__NOTE: OUT OF DATE__
Argument   |  Meaning   |   Example 
--- | --- | --- | ---
gws_path |  Path to the tabixed and bgzipped summary statistic that is going to be filtered, annotated and compared. Required argument. | path_to_summary_statistic/summary_statistic.tsv.gz 
--column-labels | Specify summary file column names. Columns specified are: (chrom, pos, ref, alt, pval, beta) . Default is '#chrom pos ref alt pval beta'. | --column-labels CHROM POS REF ALT PVAL BETA 
--sign-treshold | significance threshold for variants | --sign-treshold 5e-8 
--prefix | a prefix for all of the output files. Useful in cases where there might be confusion between processes running in the same folder. A dot is inserted after the prefix if it is passed. | --prefix diabetes_r4 
--group | supplying this flag results in the variants being grouped into groups | --group 
--grouping-method | grouping method used if --group flag is supplied. options are 'simple', i.e. grouping based on range from gws variants, or 'ld', i.e. grouping using plink --clump |  --grouping-method simple \| ld \| cred 
--locus-width-kb | group widths in kb. In case of ld clumping, the value is supplied to plink --clump-kb. | --locus-width-kb 500 
--alt-sign-treshold | optional alternate significance threshold for including less significant variants into groups. | --alt-sign-treshold 5e-6 
--ld-panel-path | path to the LD panel, without panel file suffix. LD panel must be in plink's .bed format, as a single file. Accompanying .bim and .fam files must be in the same directory. | --ld-panel-path path_to_panel/plink_file 
--ld-r2 | plink clump-r2 argument, default 0.4 | --ld-r2 0.7 
--dynamic-r2-chisq | If flag is passed, r2 threshold is set per peak so that leadvar_chisq*r2=value (default 5). | --dynamic-r2-chisq 
--ld-api | choose which LD calculation method you want to use. `online` requires no ld panel or plink usage. `plink` uses plink to calculate LD. | --ld-api plink \| online 
--plink-memory | plink --memory argument. Default 12000 | --plink-memory 16000 
--overlap | If this flag is supplied, the groups of gws variants are allowed to overlap, i.e. a single variant can appear multiple times in different groups. | --overlap 
--ignore-region| One can make the script ignore a given region in the genome, e.g. to remove HLA region from the results. The region is given in "CHR:START-END"-format. | --ignore-region 6:1-100000000 
--credible-set-file| Add SuSiE credible sets, listed in a file of .snp files. One row per .snp file.| --credible-set-file file_containing_susie_snp_files 
--previous-release-path | path to previous release summary statistic. Must be tabixed. Required for annotation. | --previous-release-path path/prev_rel_pheno.gz 
--gnomad-path | path to gnomAD v4  genome annotation file.| --gnomad-path gnomad4_path/annotation.txt.gz
--finngen-path | Path to FinnGen annotation file, containing e.g. most severe consequence and corresponding gene of the variants | --finngen-path path_to_file/annotation.tsv.gz 
--functional-path | File path to functional annotation file | --functional-path path_to_file/annotation.tsv.gz 
--use-gwascatalog | Add flag to compare results against GWAS Catalog associations | --use-gwascatalog 
--custom-dataresource | Compare against associations defined in an additional file. | --custom-dataresource file.tsv 
--report-out | comparison output file, default 'report_out.csv'. The final output of the script | --report-out output.tsv 
--gwascatalog-pval | P-value to use for filtering results from GWAS Catalog. default 5e-8 | --gwascatalog-pval 5e-6 
--gwascatalog-width-kb | Buffer outside gws variants that is searched from GWAS Catalog, in kilobases. Default 25  | --gwascatalog-width-kb 50 
--gwascatalog-threads | Number of concurrent queries to gwasgatalog API. Default 4. Increase to speed up gwascatalog comparison. | --gwascatalog-threads 8 
--gwascatalog-allele-file | compressed & tabixed VCF file which contains alleles for rsids. Required for `gwas` and `local` gwascatalog dbs. Must be matching build and version with GWAS Catalog. Currently hg38.12 b153. | --gwascatalog-allele-file
--extra-cols | Include additional columns from the summary file in the analysis, for example effect size, rsid or allele frequencies. Values are added both to variant reports and group reports, with names like lead_COLNAME in group report. | --extra-cols sebeta rsid maf_cases maf_controls 
--top-report-out | Name of per-group aggregated report output | --top-report-out top_report.csv 
--efo-traits | specific traits that you want to concentrate on the top level locus report. Other found traits will be reported on a separate column from these. Use Experimental Factor Ontology codes. | --efo-traits EFO_1 EFO_2 EFO_3 EFO_4 
--local-gwascatalog | File path to gwas catalog downloadable associations with mapped ontologies. | --local-gwascatalog gwascatalog-associations-with-ontologies.tsv 
--db | Choose which comparison database to use: GWAS Catalog proper, GWAS Catalog's summary statistic api, or a local copy of GWAS Catalog. With local copy, you need to supply the --local-gwascatalog filepath | --db gwas \| summary_stats \| local 

###  4.1.2. <a name='Example'></a>Example

For example, here's a simple bash script for running the tool with some of the more relevant parameters:
```
#! /bin/bash

file=path_to_input/summary_statistic.gz                     # input file
sig_p='5e-8'                                                # significance threshold
sig_p2='0.001'                                              # alternate significance threshold, used with grouping
prefix='analysis_prefix'                                    # output file prefix
grouping_method=cred                                        # grouping method, one of [simple, ld, cred]
loc_w=1500                                                  # range for grouping, in kilobases
ld_panel=path_to_ld/ld_panel                                # the path to the .bed LD panel, without file suffix
ld_r2=0.1                                                   # grouping (ld, cred only) LD threshold
ld_api=online                                               # Use LD calculation server to get Finnish LD. Does not require ld_panel or plink_mem!!
plink_mem=14000                                             # plink memory in MB
ignore_region='6:23000000-38000000'                         # Result region to ignore: corresponds to MHC region
credsetfile=path_to_cs/credible_set.SUSIE.snp.bgz           # SuSiE credible set output
gnomad_path=path_to_annotation/gnomad.gz                    # gnomad v4 annotation file
finngen_ann=path_to_annotation/annotation.gz                # finngen annotation file
functional_ann=path_to_annotation/functional_annotations.gz # annotation file with functional consequences
gwas_allele_path=path_to_allele_vcf/dbsnp_vcf.gz            # GWAS Catalog allele annotation file 
use_gwascatalog="--use-gwascatalog"                         # Compare against GWAS Catalog
db=gwas                                                     # GWAS Catalog database: local copy (local), normal (gwas), or summary statistic API (summ_stats)
extracols='rsid beta sebeta maf maf_cases maf_controls'     # Add additional columns to reports


python3 Scripts/main.py $file --sign-treshold $sig_p  \
        --alt-sign-treshold $sig_p2 --prefix $prefix \
        --group --grouping-method $grouping_method \
        --locus-width-kb $loc_w \
        --ld-panel-path $ld_panel --ld-r2 $ld_r2 --ld-api $ld_api \
        --plink-memory $plink_mem --ignore-region $ignore_region \
        --credible-set-file $credsetfile \
        --finngen-path $finngen_ann \
        --gnomad-path $gnomad_path
        --functional-path $functional_ann \
        --gwascatalog-allele-file  $gwas_allele_path \
        $use_gwascatalog --db $db --extra-cols $extracols
```


##  5. <a name='Outputs'></a>Outputs

__NOTE: OUT OF DATE__

### Top report
The `top_report.tsv` file contains the group-level summary of an autoreporting run. This file is mostly useful as a first step in the analysis of a  phenotype. It is a tab-separated file with one row per one credible set/group of variants. The columns are as follows: 
 
Column  |  Description  |  Example value/Formatting  
--- | --- | ---
phenotype | phenotype name | -
phenotype_abbreviation | phenotype code | -
Cases | Number of cases for this phenotype | 1234
Controls | Number of controls for this phenotype | 1234
locus_id | The locus in question, formatted from the top SNP's chromosome, position, reference and alternate alleles. In case of credible set grouping, the top SNP is the variant with the  largest PIP in that credible set. In case of LD and simple grouping, the top SNP is the variant with smallest p-value of that group/region. Most if not all release results are grouped around credible sets.| `chr1_1_C_T` for a lead variant with chromosome 1, position 1, reference allele C and alternate  allele T.
chrom | chromosome of locus | `1`
pos | lead variant position | `123456`
ref | lead variant reference allele | `A`
alt | lead variant alternate allele | `C`
pval | lead variant p-value | `5.01e-7`
start | locus start position in basepairs | `1` for a group with positions [1,2,3,4,5]
end | locus end position in basepairs | `5` for a group with positions [1,2,3,4,5]
lead_enrichment | How much the lead variant is enriched in Finnish population compared to Europe | `4.35`
lead_$COLUMN_NAME | other columns that are grabbed for the lead variant, such as allele frequencies, effect size, variant RSIDs | -
most_severe_gene | most severe gene of the lead variant | `APOE`
most_severe_consequence | most severe consequence of lead variant | `missense_variant`
gnomAD_functional_category | functional category for the variant Exome data. | `pLoF`
gnomAD_enrichment_nfsee | lead variant enrichment in Finland against NFSEE population. Exome data. | `5.1`
gnomAD_fin.AF | lead variant allele frequency in Finland. Exome data.| `0.123`
gnomAD_fin.AN | lead variant allele number in Finland. Exome data.| `123`
gnomAD_fin.AC | lead variant allele count in Finland. Exome data.| `123`
gnomAD_fin.homozygote_count | Amount of homozygote carriers in Finnish population. Exome data.| `12`
gnomAD_fet_nfsee.odds_ratio | Fischer's exact test for enrichment FIN vs NFSEE odds ratio. Exome data.| `1.15`
gnomAD_fet_nfsee.p_value | Fischer's exact test for enrichment FIN vs NFSEE p-value. Exome data.| `5.01e-3`
gnomAD_nfsee.AC | lead variant NFSEE population allele count. Exome data.| `123`
gnomAD_nfsee.AN | lead variant NFSEE population allele number. Exome data.| `123`
gnomAD_nfsee.AF | lead variant NFSEE population allele frequency. Exome data.| `0.123`
gnomAD_nfsee.homozygote_count | Amount of homozygote carriers in NFSEE population. Exome data.| `123`
cs_id | credible set id | `chr1_123456_A_C_1`
cs_size | credible set size | `5`
cs_log_bayes_factor | credible set bayes factor, log10 | `5.21`
cs_number | credible set number in its region. Numbers are not necessarily indicative of the amount of credible sets for a region. For example, a cs number of 9 does not mean there are 9 credible sets for that region. | `1`
cs_region | finemapping region | `1:1-30000001`
found_associations_strict | This column lists all of the trait associations found in GWAS Catalog for variants that are in the credible set/strict group  (strict group here means that in case of LD grouping, variants that are in higher LD than a given threshold, and have p-values lower than the significance threshold). The trait name is followed by the amount of correlation (in R²) that association had with the lead variant. If there are multiple variants associated with that trait, the largest value is chosen. | `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All of the variants associated with a trait are guaranteed to either be part of a credible set (in the case of credible set grouping), or to have LD larger than a given threshold with the top variant AND have a p-value that is significant (in case of LD grouping).
found_associations_relaxed | This column lists all of the trait associations found in GWAS Catalog for variants in the group. The trait name is followed  by the R² to lead value of the variant that had the association. If there are multiple variants associated with that trait, the largest value is chosen.  |  `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All associated variants are guaranteed to be part of this group.  
credible_set_variants | This column lists the credible set variants. The PIP and R² values are listed after the variant | `chr1_1_C_T\|0.6\|1;chr1_100_A_G\|0.2\|0.999` for variants `chr1_1_C_T` and `chr1_100_A_G`, with PIP and R² values of [0.6,0.2] and [1,0.999], respectively.
functional_variants_strict | All of the variants with a functional consequence, with the functional consequence label, gene related to that label and R² to lead variant. The variants are part of the credible set/strict group. | `chr1_1_C_T\|missense_variant\|ABC123\|0.6` for a group with one missense variant associated with gene ABC123 and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the credible set/strict group. 
functional_variants_relaxed | All of the variants with a functional consequence, with the functional consequence label, gene related to that label and R² to lead variant. | `chr1_1_C_T\|missense_variant\|ABC123\|0.6` for a group with one missense variant associated with gene ABC123 and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the group.
specific_efo_trait_associations_strict | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists only associations where the variant is in the credible set/strict group.| Same formatting as found_associations_strict 
specific_efo_trait_associations_relaxed | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists associations to all variants in the group. |  Same formatting as found_associations_relaxed
credible_set_min_r2_value | The minimum R² value to lead variant in the credible set | `0.489`


