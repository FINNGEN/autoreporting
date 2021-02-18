# Automatic hit reporting tool for FinnGen

<!-- vscode-markdown-toc -->
* 1. [Description](#Description)
* 2. [Installation](#Installation)  
    * 2.1. [Dependencies](#Dependencies)  
    * 2.2. [Installation process](#Installationprocess)  
* 3. [Resources](#Resources)  
* 4. [Usage](#Usage)  
    * 4.1. [Tool](#Completescriptmainpy)  
        * 4.1.1. [Command-line arguments](#Commandlinearguments)  
        * 4.1.2. [Example](#Example)  
    * 4.2. [Tool subroutines](#Subroutines)  
        * 4.2.1. [gws_fetch.py:](#gws_fetchpy)  
            * 4.2.1.1. [A detailed description of gws_fetch](#detailedfetch)  
        * 4.2.2. [annotate<span></span>.py](#annotatespanspanpy)  
            * 4.2.2.1. [A detailed description of annotate](#detailedannotate)  
        * 4.2.3. [compare<span></span>.py:](#comparespanspanpy)  
            * 4.2.3.1. [A detailed description of compare](#detailedcompare)  
* 5. [Outputs](#Outputs)  
* 6. [WDL pipeline](#WDLpipeline)  
    * 6.1. [WDL files](#wdlinputs)  
    * 6.2. [Running the pipeline](#wdlrun)    


<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->


#  1. <a name='Description'></a>Description

This pipeline is used to
1) Filter and group genome-wide significant variants from FinnGen summary statistics.
2) Annotate the genome-wide significant variants using gnomAD and FinnGen annotations.
3) Cross-reference the variants to previous results, e.g. GWAS Catalog or summary statistics from hand-picked studies.
<!--Currently, steps 1,3 and 4 are operational.--> 

__NOTE: currently, only files which are in build 38 are supported. This concerns all of the input files.__

__NOTE: In case you do not have access to the Finnish LD panel (requires access to refinery), use the command-line option `--ld-api online` to circumvent that and still get access to Finnish LD.__ 

#  2. <a name='Installation'></a>Installation

##  2.1. <a name='Dependencies'></a>Dependencies

Packages:
```
python 3 (version 3.6+)
pip
plink (latest binary)
ldstore (tested on 1.1)
zlib
```
Python 3 libraries:
``` 
requests
numpy
pandas
pytabix 
```

##  2.2. <a name='Installationprocess'></a>Installation process

Install dependencies (on ubuntu 18.04+):

```
sudo apt install tabix python3 python3-pip zlib1g-dev
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip 
unzip plink_linux_x86_64_latest.zip 
sudo cp plink /usr/local/bin/
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz
tar xvf ldstore_v1.1_x86_64.tgz
sudo cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/ldstore
pip3 install -r requirements.txt
```

Copy repository to folder:

```
git clone https://github.com/FinnGen/autoreporting.git
```

#  3. <a name='Resources'></a>Resources

The resources to use this tool (gnomAD & FinnGen annotations, LD panel) can be found here:
- Summary statistics: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/summary_stats/release```
- Credible sets: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/finemap```
- gnomAD genome & exome annotations: ```gs://finngen-production-library-green/autoreporting_annotations/gnomad_data/```
- FinnGen annotations: ```gs://finngen-production-library-green/finngen_R4/finngen_R4_analysis_data/annotations/```
- Functional annotations: ```gs://r4_data_west1/gnomad_functional_variants/fin_enriched_genomes_select_columns.txt.gz```
- LD panel (based on 1000 genomes data): ```gs://finngen-production-library-green/autoreporting_annotations/1kg_ld/whole_1k*```


#  4. <a name='Usage'></a>Usage

The script can be used by either calling the whole script or calling the individual scripts by themselves in the project folder.
##  4.1. <a name='Completescriptmainpy'></a>Tool

The complete autoreporting pipeline can be used by executing the ```Scripts/main.py``` script. The pipeline will then filter, group, annotate and compare the variants. Grouping method, file paths to resources, operation parameters, output file paths, comparison database and other options can be set using different command-line arguments.

```
usage: main.py [-h] [--sign-treshold SIG_TRESHOLD] [--prefix PREFIX]
               [--fetch-out FETCH_OUT] [--group]
               [--grouping-method GROUPING_METHOD]
               [--locus-width-kb LOC_WIDTH]
               [--alt-sign-treshold SIG_TRESHOLD_2]
               [--ld-panel-path LD_PANEL_PATH] [--ld-r2 LD_R2]
               [--plink-memory PLINK_MEM] [--overlap]
               [--ignore-region IGNORE_REGION]
               [--credible-set-file CRED_SET_FILE] [--ld-api LD_API_CHOICE]
               [--gnomad-genome-path GNOMAD_GENOME_PATH]
               [--gnomad-exome-path GNOMAD_EXOME_PATH] [--include-batch-freq]
               [--finngen-path FINNGEN_PATH]
               [--functional-path FUNCTIONAL_PATH]
               [--annotate-out ANNOTATE_OUT]
               [--use-gwascatalog] [--custom-dataresource CUSTOM_DATARESOURCE]
               [--check-for-ld] [--report-out REPORT_OUT]
               [--ld-report-out LD_REPORT_OUT]
               [--gwascatalog-pval GWASCATALOG_PVAL]
               [--gwascatalog-width-kb GWASCATALOG_PAD]
               [--gwascatalog-threads GWASCATALOG_THREADS]
               [--ldstore-threads LDSTORE_THREADS] [--ld-treshold LD_TRESHOLD]
               [--cache-gwas] [--column-labels CHROM POS REF ALT PVAL]
               [--extra-cols [EXTRA_COLS [EXTRA_COLS ...]]]
               [--top-report-out TOP_REPORT_OUT]
               [--strict-group-r2 STRICT_GROUP_R2]
               [--efo-codes EFO_TRAITS [EFO_TRAITS ...]]
               [--local-gwascatalog LOCALDB_PATH]
               [--db {local,gwas,summary_stats}]
               gws_fpath
```

###  4.1.1. <a name='Commandlinearguments'></a>Command-line arguments

Argument   |  Meaning   |   Example | Original script
--- | --- | --- | ---
--sign-treshold | significance threshold for variants | --sign-treshold 5e-8 | gws_fetch.py
--prefix | a prefix for all of the output and temporary files. Useful in cases where there might be confusion between processes running in the same folder. A dot is inserted after the prefix if it is passed. | --prefix diabetes_r4 | main<span></span>.py
--fetch-out | output file path for filtered and/or grouped variants. 'fetch_out.csv' by default. | --fetch-out output.tsv | gws_fetch.py
--group | supplying this flag results in the variants being grouped into groups | --group | gws_fetch.py
--grouping-method | grouping method used if --group flag is supplied. options are 'simple', i.e. grouping based on range from gws variants, or 'ld', i.e. grouping using plink --clump |  --grouping-method simple \| ld \| cred | gws_fetch.py
--locus-width-kb | group widths in kb. In case of ld clumping, the value is supplied to plink --clump-kb. | --locus-width-kb 500 | gws_fetch.py
--alt-sign-treshold | optional alternate significance threshold for including less significant variants into groups. | --alt-sign-treshold 5e-6 | gws_fetch.py
--ld-panel-path | path to the LD panel, without panel file suffix. LD panel must be in plink's .bed format, as a single file. Accompanying .bim and .fam files must be in the same directory. | --ld-panel-path path_to_panel/plink_file | gws_fetch.py
--ld-r2 | plink clump-r2 argument, default 0.4 | --ld-r2 0.7 | gws_fetch.py
--ld-api | choose which LD calculation method you want to use. `online` requires no ld panel or plink usage. `plink` uses plink to calculate LD. | --ld-api plink \| online | gws_fetch.py
--plink-memory | plink --memory argument. Default 12000 | --plink-memory 16000 | gws_fetch.py
--overlap | If this flag is supplied, the groups of gws variants are allowed to overlap, i.e. a single variant can appear multiple times in different groups. | --overlap | gws_fetch.py
--ignore-region| One can make the script ignore a given region in the genome, e.g. to remove HLA region from the results. The region is given in "CHR:START-END"-format. | --ignore-region 6:1-100000000 | gws_fetch.py
--credible-set-file| Add SuSiE credible sets, listed in a file of .snp files. One row per .snp file.| --credile-set-file file_containing_susie_snp_files | gws_fetch.py
--gnomad-genome-path | path to gnomAD genome annotation file. Must be tabixed. Required for annotation. | --gnomad-genome-path gnomad_path/gnomad_file.tsv.gz | annotate<span></span>.py
--gnomad-exome-path | path to gnomAD exome annotation file. Must be tabixed. Required for annotation. | --gnomad-exome-path gnomad_path/gnomad_file.tsv.gz | annotate<span></span>.py
--include-batch-freq | Include batch frequencies from FinnGen annotation file | --include-batch-freq | annotate<span></span>.py
--finngen-path | Path to FinnGen annotation file, containing e.g. most severe consequence and corresponding gene of the variants | --finngen-path path_to_file/annotation.tsv.gz | annotate<span></span>.py
--functional-path | File path to functional annotation file | --functional-path path_to_file/annotation.tsv.gz | annotate<span></span>.py
--annotate-out | annotation output file, default 'annotate_out.csv' | --annotate-out annotation_output.tsv | annotate<span></span>.py
--use-gwascatalog | Add flag to compare results against GWAS Catalog associations | --use-gwascatalog | compare<span></span>.py
--custom-dataresource | Compare against associations defined in an additional file. | --custom-dataresource file.tsv | compare<span></span>.py
--check-for-ld | When supplied, gws variants and summary statistics (from file or GWAS Catalog) are tested for ld using LDstore.  | --check-for-ld | compare<span></span>.py
--raport-out | comparison output file, default 'raport_out.csv'. The final output of the script, in addition to the ld_raport_out.csv, if asked for. | --raport-out raport_out.tsv | compare<span></span>.py
--ld-raport-out | ld check output file, default 'ld_raport_out.csv'. The final output of the script, in addition to the raport_out.csv. | --ld-raport-out ld_raport_out.tsv | compare<span></span>.py
--gwascatalog-pval | P-value to use for filtering results from GWAS Catalog. default 5e-8 | --gwascatalog-pval 5e-6 | compare<span></span>.py
--gwascatalog-width-kb | Buffer outside gws variants that is searched from GWAS Catalog, in kilobases. Default 25  | --gwascatalog-width-kb 50 | compare<span></span>.py
--gwascatalog-threads | Number of concurrent queries to gwasgatalog API. Default 4. Increase to speed up gwascatalog comparison. | --gwascatalog-threads 8 | compare<span></span>.py
--ldstore-threads | Number of threads to use with LDstore. At most the number of logical cores your processor has. Default 4.| --ldstore-threads 2 | compare<span></span>.py
--ld-treshold | LD threshold for LDstore, above of which summary statistic variants in ld with our variants are included. Default 0.4 | --ld-treshold 0.8 | compare<span></span>.py
--cache-gwas | Save GWAScatalog results into gwas_out_mapping.csv, from which they are read. Useful in testing. Should not be used for production runs. | --cache-gwas | compare<span></span>.py
--column-labels | Specify summary file column names. Columns specified are: (chrom, pos, ref, alt, pval) . Default is '#chrom pos ref alt pval'. | --column-labels CHROM POS REF ALT PVAL | all scripts
--extra-cols | Include additional columns from the summary file in the analysis, for example effect size, rsid or allele frequencies. Values are added both to variant reports and group reports, with names like lead_COLNAME in group report. | --extra-cols beta sebeta rsid maf_cases maf_controls | gws_fetch.py
--top-report-out | Name of per-group aggregated report output | --top-report-out top_report.csv | compare<span></span>.py
--efo-traits | specific traits that you want to concentrate on the top level locus report. Other found traits will be reported on a separate column from these. Use Experimental Factor Ontology codes. | --efo-traits EFO_1 EFO_2 EFO_3 EFO_4 | compare<span></span>.py
--local-gwascatalog | File path to gwas catalog downloadable associations with mapped ontologies. | --local-gwascatalog gwascatalog-associations-with-ontologies.tsv | compare<span></span>.py
--db | Choose which comparison database to use: GWAS Catalog proper, GWAS Catalog's summary statistic api, or a local copy of GWAS Catalog. With local copy, you need to supply the --local-gwascatalog filepath | --db gwas \| summary_stats \| local | compare<span></span>.py
gws_path |  Path to the tabixed and gzipped summary statistic that is going to be filtered, annotated and compared. Required argument. | path_to_summary_statistic/summary_statistic.tsv.gz | gws_fetch.py

The same arguments are used in the smaller scripts that the main script uses.

<!-- For example, In order to filter genome-wide significant variants, and to compare them against gwascatalog's summary statistics API, the following call can be used:
```
python3 Scripts/main.py path_to_ss/summ_stat.tsv.gz
```

Additional features, such as result grouping, can be added through the use of the aforementioned flags. -->
###  4.1.2. <a name='Example'></a>Example

For example, here's a simple bash script for running the tool with some of the more relevant parameters:
```
#! /bin/bash

file=path_to_input/summary_statistic.gz                     # input file
sig_p='5e-8'                                                # significance threshold
sig_p2='0.001'                                              # alternate significance threshold, used with grouping
prefix='analysis_prefix'                                    # output file prefix
grouping_method=ld                                          # grouping method, one of [simple, ld, cred]
loc_w=1500                                                  # range for grouping, in kilobases
ld_panel=path_to_ld/ld_panel                                # the path to the .bed LD panel, without file suffix
ld_r2=0.1                                                   # grouping (ld, cred only) LD threshold
ld_api=online                                               # Use LD calculation server to get Finnish LD. Does not require ld_panel or plink_mem!!
plink_mem=14000                                             # plink memory in MB
ignore_region='6:23000000-38000000'                         # Result region to ignore: corresponds to MHC region
credsetfile=path_to_cs/credible_set.SUSIE.snp.bgz           # SuSiE credible set output
gnomad_genome=path_to_annotation/gnomad_genomes.gz          # gnomad genome annotation file
gnomad_exome=path_to_annotation/gnomad_exomes.gz            # gnomad exome annotation file
finngen_ann=path_to_annotation/R4_annotated_variants_v1.gz  # finngen annotation file
functional_ann=path_to_annotation/functional_annotations.gz # annotation file with functional consequences
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
        --gnomad-genome-path $gnomad_genome \
        --gnomad-exome-path $gnomad_exome \
        --finngen-path $finngen_ann \
        --functional-path $functional_ann \
        $use_gwascatalog --db $db --extra-cols $extracols
```

## 4.2 <a name='Subroutines'></a>Tool subroutines

In the following subsections, the tool's main building blocks are introduced.

###  4.2.1 <a name='gws_fetchpy'></a>gws_fetch.py

```
usage: gws_fetch.py [-h] [--sign-treshold SIG_TRESHOLD] [--prefix PREFIX]
                    [--fetch-out FETCH_OUT] [--group]
                    [--grouping-method GROUPING_METHOD]
                    [--locus-width-kb LOC_WIDTH]
                    [--alt-sign-treshold SIG_TRESHOLD_2]
                    [--ld-panel-path LD_PANEL_PATH] [--ld-r2 LD_R2]
                    [--plink-memory PLINK_MEM] [--overlap]
                    [--column-labels CHROM POS REF ALT PVAL]
                    [--extra-cols [EXTRA_COLS [EXTRA_COLS ...]]]
                    [--ignore-region IGNORE_REGION]
                    [--credible-set-file CRED_SET_FILE]
                    [--ld-api LD_API_CHOICE]
                    gws_fpath
```
The gws_fetch.py script is used to filter genome-wide significant variants from the summary statistic file as well as optionally group the variants, using either location-based grouping, ld-based grouping or grouping around credible sets. The arguments used are the same as the ones in main<span></span>.py. For example, Here is a small bash script for running the gws_fetch.py script:
```
#! /bin/bash

file=path_to_input/summary_statistic.gz                     # input file
sig_p='5e-8'                                                # significance threshold
sig_p2='0.001'                                              # alternate significance threshold, used with grouping
prefix='fetch_prefix'                                       # output file prefix
grouping_method=ld                                          # grouping method, one of [simple, ld, cred]
loc_w=1500                                                  # range for grouping, in kilobases
ld_panel=path_to_ld/ld_panel                                # the path to the .bed LD panel, without file suffix
ld_api=online                                               # Use LD calculation server to get Finnish LD. Does not require ld_panel or plink_mem!!
ld_r2=0.1                                                   # grouping (ld, cred only) LD threshold
plink_mem=14000                                             # plink memory in MB
ignore_region='6:23000000-38000000'                         # Result region to ignore: corresponds to MHC region
credsetfile=path_to_cs/credible_set.SUSIE.snp.bgz           # SuSiE credible set output

python3 Scripts/gws_fetch.py $file --sign-treshold $sig_p  \
        --alt-sign-treshold $sig_p2 --prefix $prefix \
        --group --grouping-method $grouping_method \
        --locus-width-kb $loc_w \
        --ld-panel-path $ld_panel --ld-r2 $ld_r2 --ld-api $ld_api \
        --plink-memory $plink_mem --ignore-region $ignore_region \
        --credible-set-file $credsetfile
```

###  4.2.1.1. <a name='detailedfetch'></a>A detailed description of gws_fetch  
__Input__:  
gws_fpath: a FinnGen summary statistic that is gzipped and tabixed. The default column labels for the file are '#chrom','pos', 'ref', 'alt', 'pval'. These can be changed with ```--column-labels```.  
ld_panel_path (optional): a plink .bed file that is used to calculate linkage disequilibrium for the variants.  
credible_set_file: A .SUSIE.snp.bgz file, containing credible sets that were finemapped from the summary statistic.  
__Output__:  
fetch_out: a .tsv-file, with one genome-wide significant variant per row. The name of the output file can be changed using ```--fetch-out```, and a prefix can be added with ```--prefix```.  
__Filtering__:  
The tool filters the summary statistic using the p-value of the variants. The p-value threshold for filtering can be changed with ```--sign-treshold```.  
__Grouping__:
The variants can optionally be grouped into possible signals, based on one of the tree grouping methods. The grouping method is enabled with the ```--group``` flag. The grouping method can be chosen using the ```--grouping-method``` flag:
- simple (location range-based grouping): Group the variants by choosing the most significant variant and adding variants in a range around it to the group. ```--locus-width-kb``` can be used to control the range, in kilobases. ```--alt-sign-treshold``` can be used to control the p-value threshold to include variants in the groups. ```--overlap``` can be used to allow group overlap.
- ld (plink --clump): perform ld-based grouping, i.e. PLINK 1.9's --clump[1] operation on the variants. Requires the ld panel to be given to the script. ```--ld-r2``` can be used to control the plink --clump r-squared threshold. ```--plink-memory``` can be used to control the amount of memory given to plink. The same parameters as in simple grouping can be used.
- cred (credible set grouping): Assign the groups to be the credible set plus their ld-neighbors. A group always contains its own credible set, regardless of the p-value or ld threshold. In case the credible set file does not contain any credible sets, the output is an empty table. Requires the credible set file as well as the ld panel to be defined. Same parameters as in simple and ld grouping can be used.  

<!-- __Script function__:  
First, the summary statistic, a tabixed and gzipped tsv, is loaded into the script, and all variants with p-values under signifigance treshold are selected. In case no grouping is performed, they are written to a file. If grouping is performed, these variants receive an id, under column 'locus_id' in the result file, that designates which group they belong to. Grouping can be done in two ways: The variants can be simply grouped by setting a group width, or by using linkage disequilibrium to determine, how the variants are grouped. Pseudocode for the grouping by width:  
```
let group of variants G
let radius/width/distance r
while |G| > 0:
    select variant v in G with lowest p-value
    select variants v' from G, for which |v.position - v'.position| < r
    set group of v and v' to be v
    remove variants v and v' from G
```
Optionally, the groups can be set to overlap, i.e. a single variant can be included in one or more groups, if it is inside the range of the group. However, variants that are already included in a group can not form groups.  
The grouping based on linkage disequilibrium is based on plink 1.9's --clump option, which is similar in its function. Based on plink documentation, the groups are formed by taking all of the variants that are not clumped and that are inside the group's range, as well as those variants that are in ld with the group variant.[1] As such, it only differs from the simple grouping by including variants that are in ld with the group variant in the group.   -->
[1]http://zzz.bwh.harvard.edu/plink/clump.shtml

##  4.2.2. <a name='annotatespanspanpy'></a>annotate<span></span>.py

```
usage: annotate.py [-h] [--gnomad-genome-path GNOMAD_GENOME_PATH]
                   [--gnomad-exome-path GNOMAD_EXOME_PATH]
                   [--include-batch-freq] [--finngen-path FINNGEN_PATH]
                   [--functional-path FUNCTIONAL_PATH] [--prefix PREFIX]
                   [--annotate-out ANNOTATE_OUT]
                   [--column-labels CHROM POS REF ALT PVAL BETA AF AF_CASE AF_CONTROL]
                   annotate_fpath
```

The annotate<span></span>.py script is used to annotate filtered and grouped genome-wide significant variants with annotations from gnomAD and FinnGen. This script is run as a part of the main<span></span>.py, but can be run on its own as well. Here's an example bash script for running the annotate<span></span>.py script:
```

#! /bin/bash

file=fetch_output                                           # input file
prefix='annotate_prefix'                                    # output file prefix
gnomad_genome=path_to_annotation/gnomad_genomes.gz          # gnomad genome annotation file
gnomad_exome=path_to_annotation/gnomad_exomes.gz            # gnomad exome annotation file
finngen_ann=path_to_annotation/R4_annotated_variants_v1.gz  # finngen annotation file
functional_ann=path_to_annotation/functional_annotations.gz # annotation file with functional consequences

python3 Scripts/annotate.py $file --prefix $prefix \
        --gnomad-genome-path $gnomad_genome \
        --gnomad-exome-path $gnomad_exome \
        --finngen-path $finngen_ann \
        --functional-path $functional_ann \
``` 

###  4.2.2.1. <a name='detailedannotate'></a>A detailed description of annotate<span></span>.py

__input__:  
annotate_fpath: The output of gws_fetch.py  
gnomad_genome_path: A tabixed, bgzipped gnomAD genome annotation file.  
gnomad_exome_path: A tabixed, bgzipped gnomAD exome annotation file.  
finngen_path: A tabixed, bgzipped FinnGen annotation file. Batch-specific information can be included with the flag ```--include-batch-freq```. Due to file formatting changes, the FinnGen version matters. The version can be given with the ```--finngen-annotation-version``` flag. Finngen annotation files made with releases 3 or earlier should use the value 'r3', and files made with releases 4 and above should use value 'r4'. The default value is 'r3'. Given an invalid version, the script will not be able to parse the FinnGen annotation file.  
functional_path: A tabixed, bgzipped file containing the functional consequences (missense variant, pLoF etc.) for variants.  
__output__:  
annotate_out: A file with the same columns as annotate_fpath, as well as additional annotation columns from gnomAD, FinnGen and functional annotations.  

<!-- The annotate<span></span>.py-script is used to annotate the previously filtered genome-wide significant variants, using annotation files from gnomAD as well as annotation files specific to the FinnGen project. For the FinnGen annotations, specify the version/release used with the --finngen-annotation-version. Allowed values are 'r3' and 'r4'. Finngen annotations with version 3_0 or lower should be used with value 'r3', and versions 3_1 or higher (e.g. 4_1) should use the 'r4' value. The rest of the arguments are the same as in main<span></span>.py, except for annotate_fpath, which is the path to the filtered variants. For example, to annotate variants the script can be called like this:
```
python3 annotate.py variant_file_path/variants.tsv --gnomad-genome-path path_to_gnomad/gnomad_genomes.tsv.gz --gnomad-exome-path path_to_gnomad/gnomad_exomes.tsv.gz --finngen-path path_to_finngen_annotation/finngen_annotation.tsv.gz --annotate-out annotation_output.tsv
``` -->

##  4.2.3. <a name='comparespanspanpy'></a>compare<span></span>.py:

```
usage: compare.py [-h] [--sign-treshold SIG_TRESHOLD]
                  [--grouping-method GROUPING_METHOD] [--use-gwascatalog]
                  [--custom-dataresource CUSTOM_DATARESOURCE] [--check-for-ld]
                  [--plink-memory PLINK_MEM] [--ld-panel-path LD_PANEL_PATH]
                  [--prefix PREFIX] [--report-out REPORT_OUT]
                  [--ld-report-out LD_REPORT_OUT]
                  [--gwascatalog-pval GWASCATALOG_PVAL]
                  [--gwascatalog-width-kb GWASCATALOG_PAD]
                  [--gwascatalog-threads GWASCATALOG_THREADS]
                  [--ldstore-threads LDSTORE_THREADS]
                  [--ld-treshold LD_TRESHOLD] [--cache-gwas]
                  [--column-labels CHROM POS REF ALT PVAL BETA AF AF_CASE AF_CONTROL]
                  [--top-report-out TOP_REPORT_OUT]
                  [--strict-group-r2 STRICT_GROUP_R2]
                  [--efo-codes EFO_TRAITS [EFO_TRAITS ...]]
                  [--local-gwascatalog LOCALDB_PATH]
                  [--db {local,gwas,summary_stats}]
                  compare_fname
```
The compare<span></span>.py-script is used to compare the genome-wide significant variants to earlier results, either in the form of summary statistics supplied to the script or searched from GWAScatalog's summary statistic api. The optional arguments are the same as the arguments to main<span></span>.py. Here's an example bash script for running compare<span></span>.py:
```
#! /bin/bash

file=annotate_output                                        # input file
prefix='compare_prefix'                                     # output file prefix
gwascatalog_pval='5e-8'                                     # The result p-value threshold for gwas catalog hits 
gwascatalog_threads=10                                      # how many concurrent requests to gwascatalog
use_gwascatalog="--use-gwascatalog"                         # Compare results against GWAS Catalog
db=gwas                                                     # GWAS Catalog database: local copy (local), normal (gwas), or summary statistic API (summ_stats)


python3 Scripts/compare.py $file --prefix $prefix \
        --gwascatalog-pval $gwascatalog_pval \
        --gwascatalog-threads $gwascatalog_threads \
        $use_gwascatalog --db $db
```

###  4.2.3.1. <a name='detailedcompare'></a>A detailed description of compare<span></span>.py:  
__Input__:  
compare_fname: genome-wide significant variants that are filtered & grouped by gws_fetch.py and annotated by annotate<span></span>.py.  
ld_panel_path (optional): a plink .bed-file, without the suffix, that will be used by LDstore to calculate linkage disequilibrium between genome-wide significant variants and variants from other summary statistics (or GWAScatalog). Same file that's used in gws_fetch.  
summary_fpath (optional): A file containing the summary file paths. List one file path per line. Actual summary statistics must be tab-separated value files and in build 38.  
endpoint_fpath (optional):  A file containing endpoints for the summary files listed in summary_fpath. List one endpoint per line.  
custom_dataresource (optional): A tab-separated values file containing variants that you want to compare against, with one variant per row. Will be used in similar way to GWAS Catalog. Your file should have columns `chrom, pos, ref, alt, pval, beta, se, study_doi, trait, rsid, gene, notes, af, snp`. Fill missing data with `NA`.
__Output__:  
report_out: a tsv report of the variants, with each variant on its own row. If the variant has been reported in earlier studies, the phenotype and p-value for that study is announced. Variants that are novel are also reported. In case a variant associates with multiple phenotypes, all of these are reported on their own rows.  
ld_report_out: A tsv report of those variants that are in LD with external summary statistic/gwascatalog variants.  
top_report_out: A tsv report of variant groups. Information about the associated phenotypes, functional variants and credible set variants is included. More specific match information is presented in report_out. Functional variants and associations are divided into two columns, `_relaxed` and `_strict`. The `_relaxed`-columns show the functional variants/associations for the whole group, and `_strict`-columns show the functional variants/associations for 1) the credible set in case grouping around credible sets, 2) for variants in the group whose p-value is under the p-value threshold in case of other grouping methods.

__Script function__:
The comparison script takes in a filtered and annotated variant tsv file, and reports if those variants have been announced in earlier studies. The comparison can be done against GWAS Catalog, an online variant database, or against an additional dataresource. Phenotypes of interest can be added using the `--efo-codes` flag. If these phenotypes are present in the compared associations, they are presented on the top report on their own column.  
The GWAS Catalog API access can be controlled using parameters ```--gwascatalog-threads```, ```--gwascatalog-pval``` and ```--gwascatalog-width-kb```. The ```--gwascatalog-threads``` parameters decides how many concurrent connections to the GWAS Catalog are allowed. This is not limited by the number of logical cores in your computer's processor. The ```--gwascatalog-pval``` parameter sets the threshold for associations from GWAS Catalog: If the association has a smaller p-value, it is included. The ```--gwascatalog-width-kb``` parameter is not currently used.   
Optionally, genome-wide significant variants can also be tested for LD against database associations using the flag ```--check-for-ld```. Variants for which the ld value is larger than `--ld-treshold` value are reported. Same parameters that were used in gws_fetch.py, like ```--plink-memory``` and ```--ld-panel-path``` need to be used. However, these results are not incorporated to the top report.


##  5. <a name='Outputs'></a>Outputs

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
cs_number | credible set number in its region | `1`
cs_region | finemapping region | `1:1-30000001`
found_associations_strict | This column lists all of the trait associations found in GWAS Catalog for variants that are in the credible set/strict group  (strict group here means that in case of LD grouping, variants that are in higher LD than a given threshold). The trait name is followed by the amount of correlation (in R²) that association had with the lead variant. If there are multiple variants associated with that trait, the largest value is chosen. | `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All of the variants associated with a trait are guaranteed to either be part of a credible set (in the case of credible set grouping), or to have LD larger than a given threshold with the top variant (in case of LD grouping).
found_associations_relaxed | This column lists all of the trait associations found in GWAS Catalog for variants in the group. The trait name is followed  by the R² to lead value of the variant that had the association. If there are multiple variants associated with that trait, the largest value is chosen.  |  `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All associated variants are guaranteed to be part of this group.  
credible_set_variants | This column lists the credible set variants. The PIP and R² values are listed after the variant | `chr1_1_C_T\|0.6\|1;chr1_100_A_G\|0.2\|0.999` for variants `chr1_1_C_T` and `chr1_100_A_G`, with PIP and R² values of [0.6,0.2] and [1,0.999], respectively.
functional_variants_strict | All of the variants with a functional consequence, with the functional consequence label and R² to lead variant. The variants are part of the credible set/strict group. | `chr1_1_C_T\|missense_variant\|0.6` for a group with one missense variant and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the credible set/strict group. 
functional_variants_relaxed | All of the variants with a functional consequence, with the functional consequence label and R² to lead variant. The variants are part of the credible set/strict group. | `chr1_1_C_T\|missense_variant\|0.6` for a group with one missense variant and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the group.
specific_efo_trait_associations_strict | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists only associations where the variant is in the credible set/strict group.| Same formatting as found_associations_strict 
specific_efo_trait_associations_relaxed | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists associations to all variants in the group. |  Same formatting as found_associations_relaxed
credible_set_min_r2_value | The minimum R² value to lead variant in the credible set | `0.489`


##  6. <a name='WDLpipeline'></a>WDL pipeline

The WDL pipeline can be used to run a set of phenotypes as a batch job on a Cromwell server. There are two pipelines to choose from:   
- The ```autoreporting.wdl``` pipeline assigns a container for each phenotype, which creates a lot of overhead per phenotype. 
- The ```report_serial.wdl``` pipeline groups multiple phenotypes per one container, and therefore reduces the overhead per phenotype. This is mostly useful when the amount of phenotypes to process is very large, for example when running the pipeline for a whole release. This also means that running the pipeline requires more time, with the time multiplier being approximately the same as the number of phenotypes processed on one VM.


Choose the `autoreporting.wdl` pipeline for smaller runs that you need to get quickly, and `report_serial.wdl` for runs that are large and not as urgent.

### 6.1 <a name='wdlinputs'></a> WDL Files

The WDL pipeline requires an input JSON file to work correctly. Example input files can be found in the `wdl` folder. The files are named like so:  
`autoreporting_$PIPELINE_$RELEASE.json`, where the `$PIPELINE` and `$RELEASE` keywords are determined based on what pipeline and data release those files correspond to. The `serial` keyword corresponds to the `report_serial.wdl` pipeline, and the `completely_parallel` keyword corresponds to the `autoreporting.wdl` pipeline. For the releases, the only difference is in annotation resources. 

For smaller runs of the pipeline, I suggest running the `autoreporting.wdl` pipeline with a JSON file based on one of the `autoreportng_completely_parallel` files. For larger datasets, e.g over 500 phenotypes, it might be preferable to run the `report_serial.wdl` pipeline, that calculates results for multiple phenotypes on a single worker machine. For full releases, currently the only possibility is to use the `report_serial.wdl` pipeline.

The only difference between the parallel and serial JSON files is that the `autoreporting_serial_r4.json` and `autoreporting_serial_r5.json` files have an extra parameter, named `phenos_per_worker`, which determines how many phenotypes one worker machine processes. For example, if there are 1000 phenotypes and the `phenos_per_worker` parameter is set to 10, the cromwell job will spawn 100 worker machines, not 1000 worker machines like the `autoreporting.wdl` pipeline would.

### 6.2 <a name='wdlrun'></a> Running the pipeline

#### Preparing inputs
The main input for the pipeline is an input array consisting of a phenotype name, a summary statistic file path (to a google cloud bucket), and an optional SUSIE fine-mapping file. This input array is formatted as a tab-separated headerless file. You can prepare this file by hand, or if you already have a list of summary statistics and fine-mapping results, you can use the helper script `Scripts/wdl_processing_scripts/pheno_credset_array.py`. This helper script matches your summary statistic files and fine-mapping results, assuming they are named similarly.

```
usage: Create a phenotype, summary statistic file, credible set file-array from two file list files
       [-h] --phenotype-list PHENOTYPE_LIST --credset-list CREDSET_LIST
       [--phenotype-prefix PHENOTYPE_PREFIX] [--credset-prefix CREDSET_PREFIX]
       --out OUT

arguments:
  -h, --help            show this help message and exit
  --phenotype-list PHENOTYPE_LIST
                        phenotype list (e.g. output from 'gsutil ls
                        gs://pheno_folder/*.gz'
  --credset-list CREDSET_LIST
                        SuSiE credible set list (e.g. output from 'gsutil ls
                        gs://credset_folder'. NOTE: use only the SUSIE.snp.bgz-files.
  --phenotype-prefix PHENOTYPE_PREFIX
                        If the phenotypes have a prefix that is not part of
                        the phenotype name, e.g. version number use this flag
                        to include it.
  --credset-prefix CREDSET_PREFIX
                        If the credible sets have a prefix that is not part of
                        the phenotype name, e.g. version number use this flag
                        to include it.
  --out OUT             Output file, a headerless tsv file
```
#### WDL parameters

The WDL pipeline parameters are largely the same as the parameters defined for the autoreporting script. Most of these parameters can be left to be the same as in the JSON configuration files available in the folder `wdl/`, with mostly the input array and analysis-specific parameters, such as R² threshold for grouping (`autoreporting.ld_r2`) or significance threshold for including variants (`autoreporting.sign_treshold`). All of the parameters are detailed in the table below:

Parameter   |  Meaning   | Type |   Example   | Corresponding script parameter
--- | --- | --- | --- | ---
autoreporting.phenos_per_worker |  How many phenotypes to process per VM (`report_serial.wdl` only)  |   Integer   | 5 | -
autoreporting.input_array_file |  GCS (Google Cloud Storage) link to the aforementioned input array  |   String   | "gs://r5_data/autoreporting/r5_phenotype_array_2020_03_25.tsv", | -
autoreporting.docker |  URL for VM docker image  |   String   |"eu.gcr.io/finngen-refinery-dev/autorep:984100b" | -
autoreporting.memory |  VM memory in GB  |   Integer   |20 | -
autoreporting.cpus |   VM CPU count  |   Integer   |4 | -
autoreporting.gnomad_exome |   GCS link to gnomAD exome resource   |   String   |"gs://fg-datateam-analysisteam-share/gnomad/2.1/exomes/gnomad.exomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz" | `--gnomad-exome-path`
autoreporting.gnomad_genome |   GCS link to gnomAD exome resource   |   String   |"gs://fg-datateam-analysisteam-share/gnomad/2.1/genomes/gnomad.genomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz" | `--gnomad-genome-path`
autoreporting.finngen_annotation |   GCS link to FinnGen annotation resource   |   String   |"gs://r5_data/annotations/R5_annotated_variants_v1.gz" | `--finngen-path`
autoreporting.functional_annotation |   GCS link to functional annotation resource  |   String   |"gs://r4_data_west1/gnomad_functional_variants/fin_enriched_genomes_select_columns.txt.gz" | `functional-path`
autoreporting.ld_panel |   GCS link to imputation panel resource (without suffix)   |   String   |"gs://finngen-imputation-panel/sisu3/wgs_all" | `--ld-panel-path`
autoreporting.column_names |   Column names in summary statistic files. Columns are [chromosome, position, reference allele, alternate allele, p-value ]   |   Array\[String\]   |["#chrom","pos","ref","alt","pval"] | `--column-labels`
autoreporting.extra_columns |   Extra columns to include from summary statistic file than the ones presented in `column_names`, e.g. beta or MAF   |   String   |"beta maf maf_cases maf_controls rsids" | `--extra-cols`
autoreporting.efo_code_file |   GCS link to EFO to FinnGen phenotype mapping file   |   String   |"gs://r4_data_west1/autoreporting/efo_map_na.tsv" | No direct mapping, used for `--efo-codes`
autoreporting.ignore_region |   Region to ignore from the results   |   String   |"6:23000000-38000000" | `--ignore-region`
autoreporting.primary_grouping_method |   Primary grouping method, currently only makes sense to be `cred`. Primary grouping is used when the credible set file is available.   |   String   |"cred" | -
autoreporting.secondary_grouping_method |   Grouping method to use if credible set file is not available.   |   String   |"ld" | -
autoreporting.group |   Whether to group the variants or not   |   Boolean   |true | `--group`
autoreporting.overlap |   Whether variants in a group are allowed to be part of other groups   |   Boolean   | false | `--overlap`
autoreporting.include_batch_freq |   Whether to include batch-specific annotations for INFO scores   |   Boolean   |true | `--include-batch-freq`
autoreporting.sign_treshold |   Significance threshold for simple and ld grouping   |   Float   |5e-8 | `--sign-treshold`
autoreporting.alt_sign_treshold |   Alternate significance threshold for simple and ld grouping   |    Float  |1e-2 | `--alt-sign-treshold`
autoreporting.grouping_locus_width |   Grouping window in kilobases   |   Integer   |2000 | `--locus-width-kb`
autoreporting.ld_r2 |   R² threshold for LD partners   |   Float   |0.1 | `--ld-r2`
autoreporting.plink_memory |   PLINK memory in MB   |   Integer   |17000 | `--plink-memory`
autoreporting.local_gwcatalog |   GCS link to GWAS catalog resource. Needed if `db_choice` is 'local'  |   String   |"gs://r5_data/autoreporting/gwas_catalog_associations_2020-03-08.tsv" | `--local-gwascatalog`
autoreporting.db_choice |   Autoreporting GWAS catalog backend. Options: \['local','gwas','summary_stats'\]   |   String   |"local" | `--db`
autoreporting.gwascatalog_pval |   GWAS catalog p-value threshold   |   Float   |5e-8 | `--gwascatalog-pval`
autoreporting.gwascatalog_width_kb |   GWAS catalog padding in kilobases   |   Integer   | 25 | `--gwascatalog-width-pval`
autoreporting.gwascatalog_threads |   GWAS catalog API threads   |   Integer   |8 | `--gwasccatalog-threads`
autoreporting.strict_group_r2 |   R² threshold for strict grouping in case of `ld` grouping   |   Float   | 0.5 | `--strict-group-r2`
autoreporting.custom_dataresource |   GCS link to a manually curated dataresource.   |   String   |"gs://r4_data_west1/autoreporting/custom_dataresource_r4_2020_03_25.tsv" | `--custom-dataresource`

#### Running the pipeline

After the pipeline inputs are prepared, the next step is running the pipeline. This is easiest done by connecting to a cromwell server ( See instructions in https://github.com/FINNGEN/CROMWELL#access-cromwell-server-and-submit-workflows) or by using other tools for that (e.g. https://github.com/FINNGEN/CromwellInteract). In short, after the input array is created and the configuration JSON file is ready, take the WDL pipeline file and JSON file and send them to the server. The results can be acquired by copying them from the result bucket:
```
gsutil -m cp gs://[link-to-cromwell-bucket]/autoreporting/[pipeline hash]/call-report/**/*.out [destination_folder]/
``` 
<!-- 
The WDL pipeline can be used to run multiple phenotypes as a batch job on a cromwell server. This is very useful for e.g. processing all of the phenotypes in one data release. There are two different pipelines: the autoreporting.wdl is run parallel per phenotype, in that every phenotype gets its own container. This is easy to implement, but due to the large files necessary to run these pipelines (such as the LD panel), the initialization of the containers takes a substantial amount of the total processing time.  
The autoreporting_partially_serialized.wdl pipeline divides the phenotypes into smaller groups, which are then processed one group per container. This amortizes the time taken in downloading resources between the phenotypes in a group, lowering total machine time required to process all phenotypes, but possibly increasing the total time to process all phenotypes (as the time taken to process all phenotypes is the time of the slowest container). The size of these groups is changeable. The parameters in the json resource file are similarly named as the parameters in the command line.  -->


