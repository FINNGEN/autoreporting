# Automatic hit reporting tool for FINNGEN

## A tool for finding gws results from GWAS summary statistics

This pipeline is used to
1) Filter out and group gws SNPs from FINNGEN summary statistics
2) Perform finemapping on the filtered SNPs
3) Annotate the SNPs using gnoMAD and FINNGEN annotations
4) Cross-reference the SNPs to previous results, e.g. online databases and/or hand-picked results from studies

## 1: Filtering and clumping GWS variants: gws_fetch.py
The first part of the tool filters in genome-wide-significant variants (treshold configurable) and optionally groups them according to the command line settings. The grouping can be done using either grouping based on an area around top gws variants, or using PLINK's ld clumping. An alternate signifigance treshold for including variants in the groups can be specified.



## 2: Finemapping the GWS variants
TBD

## 3: Annotating the GWS variants using gnoMAD annotations and FINNGEN annotations
The third part of the tool annotates GWS variants using annotations both from gnoMAD and FINNGEN.

## 4: Comparing GWS variants to previous findings
The fourth part of the tool compares GWS variants to either external summary statistics or summary statistics in databases such as GWASCatalog.