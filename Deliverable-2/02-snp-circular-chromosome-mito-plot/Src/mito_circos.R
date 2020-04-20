#Notes: BioCircos package will not produce a graph if there is a entry that cannot be loaded


#libraries
library(readr)
library(BioCircos)
library(tidyverse)
library(data.table)
library(shiny)

# Directory setup
directory <- 'D:/data/2019-051-GamaVivian-VariantCallingAnalysis/manuscript/Deliverable-2/06-snp-circular-chromosome-mito-plot/'
setwd(directory)

#Due to differences with data and parameters used the data from GeneWiz ( DO NOT USE)
#GM01503 <- read_csv("D:/data/2019-051-GamaVivian-VariantCallingAnalysis/data/GM01503_Mitomaster query.csv", locale = locale(encoding = "ASCII"))
#GM03672 <- read_csv("D:/data/2019-051-GamaVivian-VariantCallingAnalysis/data/GM03672_Mitomaster query.csv", locale = locale(encoding = "ASCII"))
#GM13411 <- read_csv("D:/data/2019-051-GamaVivian-VariantCallingAnalysis/data/GM13411_Mitomaster query.csv", locale = locale(encoding = "ASCII"))

#Loading data and setup----------------
snpeff_gr37_snp_full <- read_csv("Input/snpeff-gr37-snp-full.csv", 
                                 col_types = cols(CHROM = col_character()))
#snpeff_gr37_indel_full <- read_csv("D:/data/2019-051-GamaVivian-VariantCallingAnalysis/data/snpeff-gr37-indel-full.csv", 
#                                   col_types = cols(CHROM = col_character()))
crossref <- read_csv("Input/crossref.csv")
color_ref <- read_csv("Input/color_ref.csv")

#removing non-MT related rows
snpeff_gr37_snp_fullmod <- snpeff_gr37_snp_full[snpeff_gr37_snp_full$CHROM == 'MT',]

#for snps mitochondria
#subset variants into samples
GM01503 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM01503_hit =='1',]

GM03672 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM03672_hit =='1',]

GM13411 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM13411_hit =='1',]

#GM01503
gm15_pos <- GM01503$`POS`
gm15_chr <- GM01503$Gene_Name
gm15_chr <- substring(gm15_chr, 4)
gm15_mut <- GM01503$Annotation_Impact
gm15 <- data.frame(chr = gm15_chr,
                   pos = gm15_pos,
                   mut = gm15_mut)

#matching the gene start and end and the length of the gene to table
idx <- match( gm15$chr, crossref$Gene)
gm15$str <- crossref$start[ idx ]
gm15$end <- crossref$end[ idx ]
gm15$chr_length <- crossref$chr_length[ idx ]

#changing the values of SNP based on the length of the gene in the genome 
#because the control region is split up (one section is at the beginning 
#of the genome and the other half is at the end), we are keeping the pos value that are under 576 nt
#and changing the values for the other SNPS
gm15$snppos <- ifelse(gm15$pos< 576, gm15$pos, gm15$chr_length-(gm15$end-gm15$pos))
gm15_2 <- subset(gm15, !is.na(str))
idx <- match( gm15_2$mut, color_ref$Mut)
gm15_2$color <- color_ref$gm15color[ idx ]
gm15_2$value <- color_ref$value[ idx ]
gm15_delet <- gm15_2[gm15_2$mut =='MODERATE',]
gm15_delet$chr <- as.character(gm15_delet$chr)

#GM03572
gm03_pos <- GM03672$`POS`
gm03_chr <- GM03672$Gene_Name
gm03_chr <- substring(gm03_chr, 4)
gm03_mut <- GM03672$Annotation_Impact
gm03 <- data.frame(chr = gm03_chr,
                   pos = gm03_pos,
                   mut = gm03_mut)

#matching the gene start and end and the length of the gene to table
idx <- match( gm03$chr, crossref$Gene)
gm03$str <- crossref$start[ idx ]
gm03$end <- crossref$end[ idx ]
gm03$chr_length <- crossref$chr_length[ idx ]

#changing the values of SNP based on the length of the gene in the genome 
#because the control region is split up (one section is at the beginning 
#of the genome and the other half is at the end), we are keeping the pos value that are under 576 nt
#and changing the values for the other 
gm03$snppos <- ifelse(gm03$pos< 576, gm03$pos, gm03$chr_length-(gm03$end-gm03$pos))
gm03_2 <- subset(gm03, !is.na(str))
idx <- match( gm03_2$mut, color_ref$Mut)
gm03_2$color <- color_ref$gm03color[ idx ]
gm03_2$value <- color_ref$value[ idx ]
gm03_delet <- gm03_2[gm03_2$mut=='MODERATE',]
gm03_delet$chr <- as.character(gm03_delet$chr)

#GM13411
gm13_pos <- GM13411$`POS`
gm13_chr <- GM13411$Gene_Name
gm13_chr <- substring(gm13_chr, 4)
gm13_mut <- GM13411$Annotation_Impact
gm13 <- data.frame(chr = gm13_chr,
                   pos = gm13_pos,
                   mut = gm13_mut)

#matching the gene start and end and the length of the gene to table
idx <- match( gm13$chr, crossref$Gene)
gm13$str <- crossref$start[ idx ]
gm13$end <- crossref$end[ idx ]
gm13$chr_length <- crossref$chr_length[ idx ]

#changing the values of SNP based on the length of the gene in the genome 
#because the control region is split up (one section is at the beginning 
#of the genome and the other half is at the end), we are keeping the pos value that are under 576 nt
#and changing the values for the other 
gm13$snppos <- ifelse(gm13$pos< 576, gm13$pos, gm13$chr_length-(gm13$end-gm13$pos))
gm13_2 <- subset(gm13, !is.na(str))
idx <- match( gm13_2$mut, color_ref$Mut)
gm13_2$color <- color_ref$gm13color[ idx ]
gm13_2$value <- color_ref$value[ idx ]
gm13_delet <- gm13_2[gm13_2$mut=='MODERATE',]
gm13_delet$chr <- as.character(gm13_delet$chr)

#BioCircos
myGenome = list("Control Region" = 1122,
                "12S" = 954,
                "16S" = 1559,
                "ND1" = 956,
                "ND2" = 1042,
                "CO1" = 1542,
                "CO2" = 684,
                "ATP8" = 207,
                "ATP6" = 681,
                "CO3" = 784,
                "ND3" = 346,
                "ND4L" = 297,
                "ND4" = 1378,
                "ND5" = 1812,
                "ND6" = 525,
                "CYB" = 1141)




#Note: Each sample track has two tracks. One track is to plot all the SNPS and the other one titled
#"Highlight" overlaps points from the first track to indicate this SNP is a deletion mutation 
#GM01503 SNP tracks 
tracklist = BioCircosSNPTrack('mySNPTrack',as.character(gm15_2$chr),gm15_2$snppos,gm15_2$value, colors = gm15_2$color, minRadius = 0.85, maxRadius = 0.95, range=c(0,12), shape = c('circle'))
tracklist = tracklist + BioCircosSNPTrack('Highlight', gm15_delet$chr, gm15_delet$snppos, gm15_delet$value, colors = gm15_delet$color, minRadius=0.85, maxRadius = 0.95, shape='circle', size = 6, range = c(0,12))

#GM03672 SNP tracks
tracklist = tracklist + BioCircosSNPTrack('mySNPTrack',as.character(gm03_2$chr),gm03_2$snppos,gm03_2$value, colors = gm03_2$color, minRadius = 0.7, maxRadius = 0.8, range=c(0,12))
tracklist = tracklist + BioCircosSNPTrack('Highlight', gm03_delet$chr, gm03_delet$snppos, gm03_delet$value, colors = gm03_delet$color, minRadius=0.7, maxRadius = 0.8, shape='circle', size = 6, range = c(0,12))

#GM13411 SNP tracks
tracklist = tracklist + BioCircosSNPTrack('mySNPTrack',as.character(gm13_2$chr),gm13_2$snppos,gm13_2$value, colors = gm13_2$color, minRadius = 0.55, maxRadius = 0.65, range=c(0,12))
tracklist = tracklist + BioCircosSNPTrack('Highlight', gm13_delet$chr, gm13_delet$snppos, gm13_delet$value, colors = gm13_delet$color, minRadius=0.55, maxRadius = 0.65, shape='circle', size = 6, range = c(0,12))


#Background tracks for SNP tracks 
tracklist = tracklist + BioCircosBackgroundTrack('GM15030',fillColors = '#EDF1EC',minRadius = 0.85, maxRadius = 0.95)

tracklist = tracklist + BioCircosBackgroundTrack('GM03672',fillColors = '#EDF1EC',minRadius = 0.7, maxRadius = 0.8)

tracklist = tracklist + BioCircosBackgroundTrack('GM13411',fillColors = '#EDF1EC',minRadius = 0.55, maxRadius = 0.65)

#tracklist = tracklist + BioCircosArcTrack('Impact',c('Control Region','Control Region','Control Region','16S','ND4'),c(150,243,510,1433,270),c(155,248,520,1440,276), colors = c('#FF0000'), minRadius = 1, maxRadius=1.14)

#BioCircos(tracklist, genome=myGenome)

#BioCircos(tracklist,genome = myGenome, genomeFillColor = "Blues",
#          genomeTicksScale = 50000000, genomeTicksDisplay = F)

#---------------------------------------------------------------------------------
#To export a SVG out BioCircos...need to output it via shiny and then use a tool to extract the SVG from the 
#webpage 
#I used a tool called SVG crowbar which is only compatiable with google chrome 

# Define shiny UI
shinyUi <- fluidPage(
  
  # Import microsoft Sans Serif font
  tags$head(
    # First load the sans serif font
    # Then apply it to the whole BioCircos HTML div id
    tags$style(HTML("
        @import url('//fonts.googleapis.com/css?family=Microsoft Sans Serif);
      ",
                    "#InstanceCircos {font-family: 'Microsoft Sans Serif', cursive;
        font-weight: 500; line-height: 1.1;
        color: #4d3a7d;}"
    ))
  ),
  
  # Set where the circos plot will be displayed
  BioCircosOutput("InstanceCircos", width = 746, height=577)
)

# Define shiny server
shinyServer <- function(input, output) {
  # Set the parameter of the circos plot
  output$InstanceCircos <- renderBioCircos({
    # Create a plot with a background track
    BioCircos(tracklist,genome = myGenome, genomeFillColor = "Blues",
              genomeTicksScale = 50000000, genomeTicksDisplay = F)
  })
}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)

