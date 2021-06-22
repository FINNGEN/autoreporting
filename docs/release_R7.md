# FINNGEN AUTOREPORTING PIPELINE RELEASE 7 NOTES  
Autoreporting pipeline, release 7
date: 2021-06-22

## File structure  
### Data  

#### Group-level reports  
Path:  
autoreporting/group_reports/
|File|Description|
|---|---| 
| PHENOTYPE.top.out | Group-level report for PHENOTYPE |

#### Variant-level reports
Path:  
autoreporting/variant_reports/

|File|Description|
|---|---| 
| PHENOTYPE.report.out | Variant-level report for PHENOTYPE |

### Documentation
|File|Description|
|---|---| 
| finngen_R7_autoreporting.pdf | This data release documentation |

## Overview

The reports were built using the SUSIE fine-mapping summaries.

Each report contains one or multiple groups, one per finemapped credible set for that phenotype. A group contains a credible set, as well as all variants that are in high LD with the credible set lead variant (LD partner). The credible set lead variant is the variant of the credible set with highest posterior inclusion probability (PIP). Depending on the quality of the credible set (good_cs true/false), the full credible set or only its lead variant is reported (as in SUSIE fine-mapping summaries). LD partners have their possible credible set status stripped, with only their credible set number kept if it was available.

Two sets of results provided for each phenotype. Variants were extracted from summary statistics using SUSIE fine-mapping credible sets. Variants in LD with these credible sets were also extracted from the summary statistics using a Finnish LD panel. The variants were then annotated and compared against GWAS Catalog. 

The tool has roughly 4 steps:

### 1. Variant filtering and grouping

Variants were restricted to fine-mapped credible sets. Variants not in credible sets were added to each locus, if 
- They were closer than 2MB of the most causal variant in the credible set, and
- They had a p-value smaller than 0.01, and
- They correlated with the group lead variant enough. The correlation threshold was calculated separately for each peak, and was set so that R^2 threshold * inverse chisquare function( lead variant p-value ) = 5.

### 2. Variant annotation

Variants were annotated with the following information:
- FinnGen R7 annotation file `gs://finngen-production-library-green/finngen_R7/finngen_R7_analysis_data/annotations/R7_annotated_variants_v1.gz`
- gnomAD exome annotation data `gs://finngen-production-library-green/autoreporting_annotations/gnomad_data/gnomad.exomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz`
- gnomAD genome annotation data `gs://finngen-production-library-green/autoreporting_annotations/gnomad_data/gnomad.genomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz`
- FinnGen R6 p-value and effect size for that variant for the same phenotype
- Coding variant information `gs://finngen-production-library-green/autoreporting_annotations/gnomad_data/fin_enriched_genomes_select_columns.txt.gz`

NOTE: In this release, variant functional consequence and related gene were taken from finngen annotation file.

### 3. Variant comparison against GWAS Catalog

Variants were searched for in GWAS Catalog. A variant was considered to match an association if
- The variant had matching chromosome sequence and chromosomal position,
- The variant was biallelic and the alleles matched,
- The p-value associated with the association phenotype was smaller than 5e-8.

### 4. Aggregation of variants into groups

Variants were aggregated per credible set.

### Output
The autoreporting tool creates two types of files: A variant report, which contains the fine-mapping results, as well as variants that were in high enough LD with the most causal variant in the credible set. These files contain one variant per row. The second type of file, a group report, aggregates these fine-mapping results, showing one credible set per row. Coding variants, found associations and credible set variants are aggregated per credible set.  
The columns for the group report can be seen in the table below.  

### Group report columns

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
found_associations_strict | This column lists all of the trait associations found in GWAS Catalog for variants that are in the credible set/strict group  (strict group here means that in case of LD grouping, variants that are in higher LD than a given threshold). The trait name is followed by the amount of correlation (in R²) that association had with the lead variant. If there are multiple variants associated with that trait, the largest value is chosen. | `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All of the variants associated with a trait are guaranteed to either be part of a credible set (in the case of credible set grouping), or to have LD larger than a given threshold with the top variant (in case of LD grouping).
found_associations_relaxed | This column lists all of the trait associations found in GWAS Catalog for variants in the group. The trait name is followed  by the R² to lead value of the variant that had the association. If there are multiple variants associated with that trait, the largest value is chosen.  |  `trait1\|1;trait2\|0.8` Given a group that has been associated with traits `trait1` and `trait2`, `trait1` association is in the top SNP and `trait2`  association with variants that have R² of [0.5,0.8] with top SNP. All associated variants are guaranteed to be part of this group.  
credible_set_variants | This column lists the credible set variants. The PIP and R² values are listed after the variant | `chr1_1_C_T\|0.6\|1;chr1_100_A_G\|0.2\|0.999` for variants `chr1_1_C_T` and `chr1_100_A_G`, with PIP and R² values of [0.6,0.2] and [1,0.999], respectively.
functional_variants_strict | All of the variants with a functional consequence, with the functional consequence label, gene related to that label and R² to lead variant. The variants are part of the credible set/strict group. | `chr1_1_C_T\|missense_variant\|ABC123\|0.6` for a group with one missense variant associated with gene ABC123 and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the credible set/strict group. 
functional_variants_relaxed | All of the variants with a functional consequence, with the functional consequence label, gene related to that label and R² to lead variant. | `chr1_1_C_T\|missense_variant\|ABC123\|0.6` for a group with one missense variant associated with gene ABC123 and R² to lead variant of `0.6`. All listed variants are guaranteed to be part of the group.
specific_efo_trait_associations_strict | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists only associations where the variant is in the credible set/strict group.| Same formatting as found_associations_strict 
specific_efo_trait_associations_relaxed | If specific traits were given to the script(e.g. equivalent EFO codes to the phenotype in question), any trait associations correspoding to those traits are listed here. This column lists associations to all variants in the group. |  Same formatting as found_associations_relaxed
credible_set_min_r2_value | The minimum R² value to lead variant in the credible set | `0.489`


### Parameters that were used:  
| parameter | value |  
|-----------|-------|
|docker image | eu.gcr.io/finngen-refinery-dev/autorep:c5c9c1d|
|gnomad exome data | gs://fg-datateam-analysisteam-share/gnomad/2.1/exomes/gnomad.exomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz |
|gnomad genome data | gs://fg-datateam-analysisteam-share/gnomad/2.1/genomes/gnomad.genomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz |
|LD panel, Finnish | gs://finngen-imputation-panel/sisu3/wgs_all.bed |
|FinnGen annotation file | gs://r7_data/annotations/R7_annotated_variants_v0.gz |
|functional annotation file | gs://finngen_commons/annotations/fin_enriched_genomes_select_columns.txt.gz |
|gwascatalog database | gs://r7_data/autoreporting/inputs/gwas-catalog-associations_ontology-annotated_2021_04_20.tsv |
|group max width (kilobases) | 2000 |
|group ld threshold |     Variable per peak (inv chisq(lead variant pval)*r2=5) |
|plink memory | 17000 |
|gwascatalog p-value threshold | 5e-8 |
|grouping | true |
|grouping method | cred |
|group overlap | false |

### Data description
Total number of phenotypes with results: 1409

### Programs used
See https://github.com/FINNGEN/autoreporting


### Formats
All results are tsv-files.
