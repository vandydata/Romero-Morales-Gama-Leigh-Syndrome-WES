# Libraries
library(splitstackshape)
library(tidyverse)
library(openxlsx)
library(biomaRt)

# Set directory
directory <-'D:/data/2019-051-GamaVivian-VariantCallingAnalysis/manuscript/Deliverable-1/02-GeneWiz_reproduce/'
setwd(directory)

# Load data
snpEff_processed_snp <- read_delim("Input/snpEff_processed.vcf", 
                                   "\t", escape_double = FALSE, col_types = cols(CHROM = col_character()), 
                                   trim_ws = TRUE)
snpEff_processed_indel<-read_delim("Input/snpEff_processed_indel.vcf",
                                   "\t", escape_double=FALSE, col_types = cols(CHROM = col_character()), trim_ws = TRUE)
ensembl_to_uniprot<-readr::read_tsv("Input/ensembl_to_uniprot")
snpeff_data<-snpEff_processed_snp
snpeff_data_indel<-snpEff_processed_indel

#changing column names to be readable
colnames(snpeff_data)[7]<-c('Allele')
colnames(snpeff_data)[8]<-c('Annotation')
colnames(snpeff_data)[9]<-c('Annotation_Impact')
colnames(snpeff_data)[10]<-c('Gene_Name')
colnames(snpeff_data)[11]<-c('Ensembl_ID')
colnames(snpeff_data)[14]<-c('GM01503_GT')
colnames(snpeff_data)[15]<-c('GM03672_GT')
colnames(snpeff_data)[16]<-c('GM13411_GT')

#removing unncessary columns
snpeff_data[12:13]<-list(NULL)

#removing any ./. genotypes from dataset, ./. means the genotype is unknown
snpeff_data$GM01503_GT[snpeff_data$GM01503_GT=='./.']<-'0/0'
snpeff_data$GM03672_GT[snpeff_data$GM03672_GT=='./.']<-'0/0'
snpeff_data$GM13411_GT[snpeff_data$GM13411_GT=='./.']<-'0/0'

#assigning hit value 
snpeff_data<-transform(snpeff_data, GM01503_hit=ifelse(GM01503_GT=='0/0',0,1))
snpeff_data<-transform(snpeff_data, GM03672_hit=ifelse(GM03672_GT=='0/0',0,1))
snpeff_data<-transform(snpeff_data, GM13411_hit=ifelse(GM13411_GT=='0/0',0,1))
snpeff_data<-transform(snpeff_data,number_of_samples_that_have_variant
                       =GM01503_hit+GM03672_hit+GM13411_hit)

#Using ZNF717 to test if we can get the same grouping counts 
snpeff_data%>% group_by(Gene_Name,number_of_samples_that_have_variant)%>% count() -> ZNF717_test
ZNF717_test[ZNF717_test$Gene_Name=='ZNF717',]

# Using Biomart to gather information on each gene 
# switching to Gr37 database
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", 'gene_biotype', 'uniprot_gn_id','description'),
                  filters = "ensembl_gene_id",
                  values = snpeff_data$Ensembl_ID,
                  mart = ensembl )



#biomaRt related data to put into the table
idx <- match( snpeff_data$Ensembl_ID, genemap$ensembl_gene_id )
snpeff_data$entrezgene <- genemap$entrezgene[ idx ]
snpeff_data$Gene_description<-genemap$description[idx]

#ensembl to uniprot data (biomart by different author)
idx<-match(snpeff_data$Ensembl_ID,ensembl_to_uniprot$`Ensembl gene ID`)
snpeff_data$UniProt_ID<-ensembl_to_uniprot$`UniProt ID(supplied by UniProt)`[idx]

#splitting the gene description to obtain the gene name(not symbol)
snpeff_data$Gene_description<-sapply(strsplit(snpeff_data$Gene_description,'[',fixed=TRUE),'[',1)
write.csv(snpeff_data,'Output/snpeff-gr37-snp-full.csv')

#hyperlinks
url<-paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",snpeff_data$Gene_Name)
url<-paste0('=HYPERLINK("', url, '"',',','"',snpeff_data$Gene_Name,'")')
snpeff_data$url<-url

#rearranging columns order to similar to GeneWiz report and changing the names
snpeff_data$Gene_Name<-NULL
colnames(snpeff_data)[colnames(snpeff_data)=="url"]<-"Gene_Symbol"
colnames(snpeff_data)[colnames(snpeff_data)=='Gene_description']<-"Gene_Name"
snpeff_firstpass<-snpeff_data[,c(21,18,19,20,17,14:16,1:13)]
write.xlsx(snpeff_firstpass,"Output/snpeff_snp_Gr37.xlsx")

#GeneWiz reproduce report for Indels
genewiz_rereport_snp<-snpeff_firstpass%>%group_by(Gene_Symbol,entrezgene,Gene_Name,
                                                  UniProt_ID,number_of_samples_that_have_variant,Ensembl_ID) %>% count()
#write genewiz reproduced report with additional items
write.xlsx(genewiz_rereport_snp,"Output/genewiz_report_snpeff_snp_Gr37.xlsx")


###-----Procedure for indel----- 

#changing column names to be readable
colnames(snpeff_data_indel)[7]<-c('Allele')
colnames(snpeff_data_indel)[8]<-c('Annotation')
colnames(snpeff_data_indel)[9]<-c('Annotation_Impact')
colnames(snpeff_data_indel)[10]<-c('Gene_Name')
colnames(snpeff_data_indel)[11]<-c('Ensembl_ID')
colnames(snpeff_data_indel)[14]<-c('GM01503_GT')
colnames(snpeff_data_indel)[15]<-c('GM03672_GT')
colnames(snpeff_data_indel)[16]<-c('GM13411_GT')

#removing unncessary columns
snpeff_data_indel[12:13]<-list(NULL)

#removing any ./. genotypes from dataset, ./. means the genotype is unknown
snpeff_data_indel$GM01503_GT[snpeff_data_indel$GM01503_GT=='./.']<-'0/0'
snpeff_data_indel$GM03672_GT[snpeff_data_indel$GM03672_GT=='./.']<-'0/0'
snpeff_data_indel$GM13411_GT[snpeff_data_indel$GM13411_GT=='./.']<-'0/0'

#assigning hit value 
snpeff_data_indel<-transform(snpeff_data_indel, GM01503_hit=ifelse(GM01503_GT=='0/0',0,1))
snpeff_data_indel<-transform(snpeff_data_indel, GM03672_hit=ifelse(GM03672_GT=='0/0',0,1))
snpeff_data_indel<-transform(snpeff_data_indel, GM13411_hit=ifelse(GM13411_GT=='0/0',0,1))
snpeff_data_indel<-transform(snpeff_data_indel,number_of_samples_that_have_variant
                             =GM01503_hit+GM03672_hit+GM13411_hit)

#Using ZNF717 to test if we can get the same grouping counts 
snpeff_data_indel%>% group_by(Gene_Name,number_of_samples_that_have_variant)%>% count() -> ZNF717_test
ZNF717_test[ZNF717_test$Gene_Name=='ZNF717',]


# #switching to Gr37 database
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", 'gene_biotype', 'uniprot_gn_id','description'),
                  filters = "ensembl_gene_id",
                  values = snpeff_data_indel$Ensembl_ID,
                  mart = ensembl )

#biomaRt related data to put into the table
idx <- match( snpeff_data_indel$Ensembl_ID, genemap$ensembl_gene_id )
snpeff_data_indel$entrezgene <- genemap$entrezgene[ idx ]
snpeff_data_indel$Gene_description<-genemap$description[idx]

#ensembl to uniprot data (biomart by different author)
idx<-match(snpeff_data_indel$Ensembl_ID,ensembl_to_uniprot$`Ensembl gene ID`)
snpeff_data_indel$UniProt_ID<-ensembl_to_uniprot$`UniProt ID(supplied by UniProt)`[idx]

#splitting the gene description to obtain the gene name(not symbol)
snpeff_data_indel$Gene_description<-sapply(strsplit(snpeff_data_indel$Gene_description,'[',fixed=TRUE),'[',1)
write.csv(snpeff_data_indel,'Output/snpeff-gr37-indel-full.csv')

#hyperlinks
url<-paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",snpeff_data_indel$Gene_Name)
url<-paste0('=HYPERLINK("', url, '"',',','"',snpeff_data_indel$Gene_Name,'")')
snpeff_data_indel$url<-url

#rearranging columns order to similar to GeneWiz report and changing the names
snpeff_data_indel$Gene_Name<-NULL
colnames(snpeff_data_indel)[colnames(snpeff_data_indel)=="url"]<-"Gene_Symbol"
colnames(snpeff_data_indel)[colnames(snpeff_data_indel)=='Gene_description']<-"Gene_Name"
snpeff_indel_firstpass<-snpeff_data_indel[,c(21,18,19,20,17,14:16,1:13)]
write.xlsx(snpeff_indel_firstpass,"Output/snpeff_indel_Gr37.xlsx")

#GeneWiz reproduce report for Indels
genewiz_rereport_indel<-snpeff_indel_firstpass%>%group_by(Gene_Symbol,entrezgene,Gene_Name,
                                                          UniProt_ID,number_of_samples_that_have_variant,Ensembl_ID) %>% count()
#write genewiz reproduced report with additional items
write.xlsx(genewiz_rereport_indel,"Output/genewiz_report_snpeff_indel_Gr37.xlsx")

