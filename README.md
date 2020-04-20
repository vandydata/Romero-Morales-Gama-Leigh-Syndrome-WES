# Modeling the effect of Leigh syndrome-associated mutations using human cerebral organoids

Authors: Alejandra I. Romero-Morales, Anuj Rastogi, Hoor Temuri, Megan Rasmussen, Gregory Scott McElroy, Lawrence Hsu, Paula M. Almonacid, Bryan Millis, Navdeep S Chandel, Jean-Philippe Cartailler and Vivian Gama

Please cite: TBD

## Data and metadata

All raw data has been deposited via GEO:

- SRA accession: 

  [PRJNA626388]: https://www.ncbi.nlm.nih.gov/sra/PRJNA626388

   (A reviewer link is available upon request)

- Release date: 2020-08-01 (tentative until manuscript is accepted for publication)

## Source code

### Variant analysis (deliverable 1)

####  01-SnpEff:

##### Input:

SNP: `Deliverable-1/01-Preprocessing/Input/all.samples.filtered.snp.nonsynonymous.ann.vcf`

Indel: `Deliverable-1/01-Preprocesing/Input/all.samples.filtered.indel.nonsynonymous.ann.vcf`

##### Procedure:

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
`Deliverable-1/02-GeneWiz_reproduce/Input/snpEff_processed_indel.vcf`,
`Deliverable-1/02-GeneWiz_reproduce/Input/ensembl_to_uniprot`

##### Procedure:

Run `Deliverable-1/02-GeneWiz_reproduce/Src/snpeff_reproduce_GeneWiz_Variant_Hit.R`

##### Output:

Location: `Deliverable-1/02-GeneWiz_reproduce/Output/`

Items:

`genewiz_report_snpeff_indel_Gr37.xlsx` and ` genewiz_report_snpeff_snp_Gr37.xlsx` - reproduction of the GeneWiz report from raw data

`snpeff_indel_Gr37.xlsx` and `snpeff_snp_Gr37.xlsx` - additional information CDS added on top of the genewiz report

`snpeff-gr37-indel-full.csv` and `snpeff-gr37-snp-full.csv`- Input data generated for other deliverable.

### Figure 1B
TBD

### Figure 1C
TBD

### Figure 1D
TBD


