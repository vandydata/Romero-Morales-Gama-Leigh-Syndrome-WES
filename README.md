# Modeling the effect of Leigh syndrome-associated mutations using human cerebral organoids

Authors: Alejandra I. Romero-Morales, Anuj Rastogi, Hoor Temuri, Megan Rasmussen, Gregory Scott McElroy, Lawrence Hsu, Paula M. Almonacid, Bryan Millis, Navdeep S Chandel, Jean-Philippe Cartailler and Vivian Gama

Please cite: TBD

## Data and metadata

All raw data has been deposited via GEO:

- SRA accession: [PRJNA626388]( https://www.ncbi.nlm.nih.gov/sra/PRJNA626388)  (A reviewer link is available upon request)
- Release date: 2020-08-01 (tentative until manuscript is accepted for publication)

## Source code

### Variant analysis (deliverable 1)

####  01-SnpEff:

##### Input:

SNP: `Deliverable-1/01-Preprocessing/Input/all.samples.filtered.snp.nonsynonymous.ann.vcf`

Indel: `Deliverable-1/01-Preprocesing/Input/all.samples.filtered.indel.nonsynonymous.ann.vcf`

##### Code:

Used snpSift to extract the columns in both the SNP and indel VCF files before downstream analysis:

SNP file:

```
java -jar SnpSift.jar extractFields all.samples.filtered.snp.nonsynonymous.ann.vcf CHROM POS ID REF ALT QUAL "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "GEN[GM01503].GT" "GEN[GM03672].GT" "GEN[GM13411].GT" > snpEff_processed.vcf
```

Indel file:

```
java -jar SnpSift.jar extractFields all.samples.filtered.indel.nonsynonymous.ann.vcf CHROM POS ID REF ALT QUAL "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "GEN[GM01503].GT" "GEN[GM03672].GT" "GEN[GM13411].GT" > snpEff_processed_indel.vcf
```

##### Output:

SNP: `Deliverable-1/01-Preprocessing/Output/snpEff_processed.vcf`
Indel: `Deliverable-1/01-Preprocessing/Output/snpEff_processed_indel.vcf`

#### 02-GeneWiz_reproduce

##### Input:

`Deliverable-1/02-GeneWiz_reproduce/Input/snpEff_processed.vcf`
`Deliverable-1/02-GeneWiz_reproduce/Input/snpEff_processed_indel.vcf`
`Deliverable-1/02-GeneWiz_reproduce/Input/ensembl_to_uniprot`

##### Code:

Run `Deliverable-1/02-GeneWiz_reproduce/Src/snpeff_reproduce_GeneWiz_Variant_Hit.R`

##### Output:

Location: `Deliverable-1/02-GeneWiz_reproduce/Output/`

Items:

`genewiz_report_snpeff_indel_Gr37.xlsx` and ` genewiz_report_snpeff_snp_Gr37.xlsx` - reproduction of the GeneWiz report from raw data

`snpeff_indel_Gr37.xlsx` and `snpeff_snp_Gr37.xlsx` - additional information CDS added on top of the genewiz report

`snpeff-gr37-indel-full.csv` and `snpeff-gr37-snp-full.csv`- Input data generated for other deliverable.

### Figure 1B (Indel circular plot)

#### Input: 

`Deliverable-2/01-snp-indel-circular-plot/Input/snpeff-gr37-indel-full.csv` 

#### Code:

1. Run `Deliverable-2/01-snp-indel-circular-plot/Src/biocircos.R`
2. Save plots.

Note -  three genes manually using Inkscape 0.92. Ticks have been drawn in.

| Gene Name  | chrom | end       | start     |
| ---------- | ----- | --------- | --------- |
| AC074212.3 | 19    | 46236509  | 46267792  |
| AC012123.1 | 18    | 30354376  | 30349758  |
| AC096644.1 | 1     | 220608023 | 220603286 |

#### Output:

Location: `Deliverable-2/01-snp-indel-circular-plot/Output/`

SVG formatted circular plot:`indel_circosplot_final.svg`

### Figure 1C  (SNP circular plot)

#### Input: 

`Deliverable-2/01-snp-indel-circular-plot/Input/snpeff-gr37-snp-full.csv`

#### Code:

1. Run `Deliverable-2/01-snp-indel-circular-plot/Src/biocircos.R`
2. Save plots.

#### Output:

Location: `Deliverable-2/01-snp-indel-circular-plot/Output/`

SVG formatted circular plot: `snp_circosplot_final.svg`


### Figure 1D (Mitochondrial circular plot)

#### Input:

`Deliverable-2/02-snp-circular-chromosome-mito-plot/Input/snpeff-gr37-snp-full.csv`

`Deliverable-2/02-snp-circular-chromosome-mito-plot/Input/color_ref.csv`

`Deliverable-2/02-snp-circular-chromosome-mito-plot/Input/crossref.csv`

#### Code:

1. Run `Deliverable-2/02-snp-circular-chromosome-mito-plot/Src/mito_circos.R`
2. Save plots

#### Output:

Location: `Deliverable-2/02-snp-circular-chromosome-mito-plot/Output/`

SVG formatted circular plot: `mito_circoplot.svg`

