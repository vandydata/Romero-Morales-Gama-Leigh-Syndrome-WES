#libraries---------------
library(BioCircos)
library(tidyverse)
library(readr)
library(biomaRt)
library(scales)
library(shiny)

# Set directory
directory <- 'D:/data/2019-051-GamaVivian-VariantCallingAnalysis/manuscript/Deliverable-2/01-snp-indel-circular-plot/'
setwd(directory) 

#Loading data and setup----------------
snpeff_gr37_snp_full <- read_csv("Input/snpeff-gr37-snp-full.csv", 
                                 col_types = cols(CHROM = col_character()))
snpeff_gr37_indel_full <- read_csv("Input/snpeff-gr37-indel-full.csv", 
                                 col_types = cols(CHROM = col_character()))

#removing MT related rows
snpeff_gr37_snp_fullmod <- snpeff_gr37_snp_full[snpeff_gr37_snp_full$CHROM != 'MT',]

snpeff_gr37_indel_fullmod <- snpeff_gr37_indel_full[snpeff_gr37_indel_full$CHROM != 'MT',]

#subset variants into samples 
GM01503 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM01503_hit =='1',]

GM03672 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM03672_hit =='1',]

GM13411 <- snpeff_gr37_snp_fullmod[snpeff_gr37_snp_fullmod$GM13411_hit =='1',]



#cross reference table to convert high qual to low
before <- c(1,2,3,4,5,6,7,8,9,10,11,12)
after <- c(12,11,10,9,8,7,6,5,4,3,2,1)
crossref <- data.frame(before,after)

#Extracting each sample variants chromosome, position and quality score
#rescale the quality scores to ensure data fits within the parameters of the package
#flip flop the values of the quality scores to make the snp tracks go outward instead of  inward
gm15_chr <- GM01503$CHROM
gm15_pos <- GM01503$POS
gm15_value <- GM01503$QUAL
gm15_value_sd <- round(rescale(gm15_value, c(1,12)),0)
idx <- match(gm15_value_sd, crossref$before)
gm15_newvalue <- crossref$after[ idx ]

gm03_chr <- GM03672$CHROM
gm03_pos <- GM03672$POS
gm03_value <- GM03672$QUAL
gm03_value_sd <- round(rescale(gm03_value, c(1,12)),0)
idx <- match(gm03_value_sd, crossref$before)
gm03_newvalue <- crossref$after[ idx ]

gm13_chr <- GM13411$CHROM
gm13_pos <- GM13411$POS
gm13_value <- GM13411$QUAL
gm13_value_sd <- round(rescale(gm13_value, c(1,12)),0)
idx <- match(gm13_value_sd, crossref$before)
gm13_newvalue <- crossref$after[ idx ]



#Map the genes start and end positions to our dataset
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", 'gene_biotype', 'uniprot_gn_id','description','start_position','end_position'),
                  filters = "ensembl_gene_id",
                  values = snpeff_gr37_snp_fullmod$Ensembl_ID,
                  mart = ensembl )

#subset the top 20 genes to map
top20 <- snpeff_gr37_snp_fullmod %>% group_by(Gene_Name,number_of_samples_that_have_variant) %>% count() 
top20 <- top20[ top20$number_of_samples_that_have_variant == 3 ,] %>% arrange(desc(n))
top20 <- top20[1:20 ,]

#map the start and end positions of top 20 genes 
idx <- match(top20$Gene_Name, genemap$hgnc_symbol)
top20$start <- genemap$start_position[idx]
top20$end <- genemap$end_position[idx]
idx <- match(top20$Gene_Name, snpeff_gr37_snp_fullmod$Gene_Name)
top20$chrom <- snpeff_gr37_snp_fullmod$CHROM[idx]

top20$length <- top20$end-top20$start
top20$chrom <- as.double(top20$chrom)
top20 <- top20 %>% arrange(desc(chrom))


#----------------------------------
#Creating the tracks 


# Arc track values used to show where the genes are 
chromosomes <- top20$chrom
start <- top20$start-1000000
end <- top20$end+1000000

#color values for snps
#create crossreference to label certain values with a color
#light color = low quality
#dark color = high quality

#GM01503 (PDH) red color
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#FF0000","#FF0000","#FF0000",'#DF0101','#DF0101','#DF0101','#B40404','#B40404','#B40404','#8A0808','#8A0808','#8A0808')
colorxross <- data.frame(value,color)
idx <- match(gm15_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm15 <- as.character(colorvalue)

#GM03672 (DLD) green color 
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#00FF00","#00FF00","#00FF00",'#01DF01','#01DF01','#01DF01','#04B404','#04B404','#04B404','#088A08','#088A08','#088A08')
colorxross <- data.frame(value,color)
idx <- match(gm03_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm03 <- as.character(colorvalue)

#GM13411 (MTMT-ATP6) orange color 
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#FF8000","#FF8000","#FF8000",'#B45F04','#B45F04','#B45F04','#8A4B08','#8A4B08','#8A4B08','#2A1B0A','#2A1B0A','#2A1B0A')
colorxross <- data.frame(value,color)
idx <- match(gm13_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm13 <- as.character(colorvalue)


#Set the snps for each sample 
tracks = BioCircosSNPTrack("testSNP1", gm15_chr, gm15_pos, gm15_newvalue, colors = colorgm15, 
                           minRadius=0.8,maxRadius = 0.9, range = c(2,12))
tracks = tracks +  BioCircosSNPTrack("testSNP1", gm03_chr, gm03_pos,  gm03_newvalue, colors = colorgm03, 
                                    minRadius=0.6,maxRadius = 0.75, range = c(2,12))
tracks = tracks + BioCircosSNPTrack("testSNP1", gm13_chr , gm13_pos, gm13_newvalue, colors = colorgm13, 
                                    minRadius=0.45,maxRadius = 0.55, range = c(2,12))

#set background tracks for the snps 
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.8, maxRadius = 0.92,
                                     borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.6, maxRadius = 0.77,
                                           borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.45, maxRadius = 0.58,
                                           borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")  

#Setting the Arc Track---------------------------------------
tracks = tracks + BioCircosArcTrack("fredTestArc", chromosomes, 
                                    starts = start, ends = end, 
                                    maxRadius = 1.15, minRadius = 1, color='#000000')

#Labeling the top 20 genes 
tracks = tracks + BioCircosTextTrack('myTextTrack', 'HNRNPCL1', size = "0.6em", opacity = 0.8, 
                               x = -.1, y =-1.25)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'IGSF3', size = "0.6em", opacity = 0.8, 
                                    x = .3, y =-1.25)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'ANKRD36C', size = "0.6em", opacity = 0.8, 
                                    x = .8, y =-0.95)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'FRG2C, ZNF717', size = "0.6em", opacity = 0.8,
                                    x = 1.1, y =-0.5)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'HLA-DRB1', size = "0.6em", opacity = 0.8, 
                                    x = 1, y =0.7)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'PRSS3', size = "0.6em", opacity = 0.8, 
                                    x = -.1, y =1.25)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'PRSS1, TRBV10-1', size = "0.6em", opacity = 0.8,
                                    x = .3, y =1.23)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'OR8U1', size = "0.6em", opacity = 0.8, 
                                    x = -.95, y =1)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'MUC6', size = "0.6em", opacity = 0.8, 
                                    x = -.8, y =1.1)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'TAS2R31, TAS2R46', size = "0.6em", opacity = 0.8, 
                                    x = -1.3, y =0.85)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'CTDSP2, KRT18', size = "0.6em", opacity = 0.8,
                                    x = -1.43, y =0.7)
tracks = tracks + BioCircosTextTrack('myTextTrack','PABPC3', size = "0.6em", opacity = 0.8, 
                                    x = -1.39, y =0.5)
tracks = tracks + BioCircosTextTrack('myTextTrack','CDC27', size = "0.6em", opacity = 0.8, 
                                    x = -1.39, y =-.55)
tracks = tracks + BioCircosTextTrack('myTextTrack','HYDIN', size = "0.6em", opacity = 0.8, 
                                    x = -1.39, y =-.4)
tracks = tracks + BioCircosTextTrack('myTextTrack','MUC16', size = "0.6em", opacity = 0.8, 
                                    x = -1.19, y =-.8)
tracks = tracks + BioCircosTextTrack('myTextTrack','FAM182B', size = "0.6em", opacity = 0.8, 
                                    x = -1.0, y =-1)


#--------------------------------------------------------------------

#To export a SVG out BioCircos...need output it via shiny and then use a tool to extract the SVG from the 
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
    BioCircos(tracks, genomeFillColor = "Spectral", yChr = T, chrPad = .012, displayGenomeBorder = F,
              genomeTicksLen = 3, genomeTicksTextSize = 0, genomeTicksScale = 5000000000,
              genomeLabelTextSize = 12, genomeLabelDy = -47, genomeTicksDisplay = F)
  })
}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)


####Indel graphic-----

# Note: There are three empty slots in the graphic. I used inkscape to manually add those into the image.
# The symbols are AC074212.3, AC012123.1, and AC096644.1

# Subset variants into samples 
GM01503 <- snpeff_gr37_indel_fullmod[snpeff_gr37_indel_fullmod$GM01503_hit =='1',]

GM03672 <- snpeff_gr37_indel_fullmod[snpeff_gr37_indel_fullmod$GM03672_hit =='1',]

GM13411 <- snpeff_gr37_indel_fullmod[snpeff_gr37_indel_fullmod$GM13411_hit =='1',]



# Cross reference table to convert high qual to low
before <- c(1,2,3,4,5,6,7,8,9,10,11,12)
after <- c(12,11,10,9,8,7,6,5,4,3,2,1)
crossref <- data.frame(before,after)

# Extracting each sample variants chromosome, position and quality score
# rescale the quality scores to ensure data fits within the parameters of the package
# flip flop the values of the quality scores to make the snp tracks go outward instead of  inward
gm15_chr <- GM01503$CHROM
gm15_pos <- GM01503$POS
gm15_value <- GM01503$QUAL
gm15_value_sd <- round(rescale(gm15_value, c(1,12)),0)
idx <- match(gm15_value_sd, crossref$before)
gm15_newvalue <- crossref$after[ idx ]

gm03_chr <- GM03672$CHROM
gm03_pos <- GM03672$POS
gm03_value <- GM03672$QUAL
gm03_value_sd <- round(rescale(gm03_value, c(1,12)),0)
idx <- match(gm03_value_sd, crossref$before)
gm03_newvalue <- crossref$after[ idx ]

gm13_chr <- GM13411$CHROM
gm13_pos <- GM13411$POS
gm13_value <- GM13411$QUAL
gm13_value_sd <- round(rescale(gm13_value, c(1,12)),0)
idx <- match(gm13_value_sd, crossref$before)
gm13_newvalue <- crossref$after[ idx ]



# Map the genes start and end positions to our dataset
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", 'gene_biotype', 'uniprot_gn_id','description','start_position','end_position'),
                  filters = "ensembl_gene_id",
                  values = snpeff_gr37_indel_fullmod$Ensembl_ID,
                  mart = ensembl )

# Subset the top 20 genes to map
top20 <- snpeff_gr37_indel_fullmod %>% group_by(Gene_Name,number_of_samples_that_have_variant) %>% count() 
top20 <- top20[ top20$number_of_samples_that_have_variant == 3 ,] %>% arrange(desc(n))
top20 <- top20[1:20 ,]

idx <- match(top20$Gene_Name, snpeff_gr37_indel_full$Gene_Name)
top20$ensemblid <- snpeff_gr37_indel_full$Ensembl_ID[idx]

# Map the start and end positions of top 20 genes 
idx <- match(top20$ensemblid, genemap$ensembl_gene_id)
top20$start <- genemap$start_position[idx]
top20$end <- genemap$end_position[idx]
idx <- match(top20$Gene_Name, snpeff_gr37_indel_fullmod$Gene_Name)
top20$chrom <- snpeff_gr37_indel_fullmod$CHROM[idx]

top20$length <- top20$end-top20$start
top20$chrom <- as.double(top20$chrom)
top20 <- top20 %>% arrange(desc(chrom))


#----------------------------------
# Creating the tracks 


# Arc track values used to show where the genes are 
chromosomes <- top20$chrom
start <- top20$start-1000000
end <- top20$end+1000000

# Color values for snps
# created a cross reference to label certain values with a color
# light color = low quality
# dark color = high quality

# GM01503 (PDH) red color
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#FF0000","#FF0000","#FF0000",'#DF0101','#DF0101','#DF0101','#B40404','#B40404','#B40404','#8A0808','#8A0808','#8A0808')
colorxross <- data.frame(value,color)
idx <- match(gm15_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm15 <- as.character(colorvalue)

# GM03672 (DLD) green color 
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#00FF00","#00FF00","#00FF00",'#01DF01','#01DF01','#01DF01','#04B404','#04B404','#04B404','#088A08','#088A08','#088A08')
colorxross <- data.frame(value,color)
idx <- match(gm03_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm03 <- as.character(colorvalue)

# GM13411 (MTMT-ATP6) orange color 
value <- c(12,11,10,9,8,7,6,5,4,3,2,1)
color <- c("#FF8000","#FF8000","#FF8000",'#B45F04','#B45F04','#B45F04','#8A4B08','#8A4B08','#8A4B08','#2A1B0A','#2A1B0A','#2A1B0A')
colorxross <- data.frame(value,color)
idx <- match(gm13_newvalue, colorxross$value)
colorvalue <- colorxross$color[idx]
colorgm13 <- as.character(colorvalue)


# Set the snps for each sample 
tracks = BioCircosSNPTrack("testSNP1", gm15_chr, gm15_pos, gm15_newvalue, colors = colorgm15, 
                           minRadius=0.8,maxRadius = 0.9, range = c(2,12))
tracks = tracks +  BioCircosSNPTrack("testSNP1", gm03_chr, gm03_pos,  gm03_newvalue, colors = colorgm03, 
                                     minRadius=0.6,maxRadius = 0.75, range = c(2,12))
tracks = tracks + BioCircosSNPTrack("testSNP1", gm13_chr , gm13_pos, gm13_newvalue, colors = colorgm13, 
                                    minRadius=0.45,maxRadius = 0.55, range = c(2,12))

# Set background tracks for the snps 
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.8, maxRadius = 0.92,
                                           borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.6, maxRadius = 0.77,
                                           borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")
tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.45, maxRadius = 0.58,
                                           borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#D8D8D8")  


# Setting the Arc Track---------------------------------------
tracks = tracks + BioCircosArcTrack("fredTestArc", chromosomes, 
                                    starts = start, ends = end, 
                                    maxRadius = 1.15, minRadius = 1, color='#000000')

# Labeling the top 20 genes 
tracks = tracks + BioCircosTextTrack('myTextTrack', 'UBXN11', size = "0.6em", opacity = 0.8, 
                                    x = .00, y =-1.25)
tracks = tracks + BioCircosTextTrack('myTextTrack', 'ALMS1', size = "0.6em", opacity = 0.8, 
                                    x = .75, y =-1)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'ANKRD36C', size = "0.6em", opacity = 0.8, 
                                    x = .8, y =-0.95)

tracks = tracks + BioCircosTextTrack('myTextTrack','ZNF717', size = "0.6em", opacity = 0.8,
                                    x = 1.1, y =-0.5)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'DSPP', size = "0.6em", opacity = 0.8, 
                                    x = 1.2, y =-0.01)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'TRBV10-1', size = "0.6em", opacity = 0.8, 
                                    x = .45, y =1.15)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'PRSS3,ADAMTSL1,ARHGAP39', size = "0.6em", opacity = 0.8,
                                    x = -.2, y =1.23)


tracks = tracks + BioCircosTextTrack('myTextTrack', 'MUC6', size = "0.6em", opacity = 0.8, 
                                    x = -0.67, y =1.2)

tracks = tracks + BioCircosTextTrack('myTextTrack', 'TAS2R43,ATN1', size = "0.6em", opacity = 0.8,
                                    x = -1.2, y =0.9)


tracks = tracks + BioCircosTextTrack('myTextTrack','PABPC3', size = "0.6em", opacity = 0.8, 
                                    x = -1.3, y =.6)

tracks = tracks + BioCircosTextTrack('myTextTrack','IGHJ6', size = "0.6em", opacity = 0.8, 
                                    x = -1.39, y =0.12)

tracks = tracks + BioCircosTextTrack('myTextTrack','HYDIN', size = "0.6em", opacity = 0.8, 
                                    x = -1.35, y =-.3)

tracks = tracks + BioCircosTextTrack('myTextTrack','CDC27', size = "0.6em", opacity = 0.8, 
                                    x = -1.3, y =-.5)

tracks = tracks + BioCircosTextTrack('myTextTrack','MUC16', size = "0.6em", opacity = 0.8, 
                                    x = -1.15, y =-.8)

shinyApp(ui = shinyUi, server = shinyServer)


#--------------------------------------------------------------------

#To export a SVG out BioCircos...need output it via shiny and then use a tool to extract the SVG from the 
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
    BioCircos(tracks, genomeFillColor = "Spectral", yChr = T, chrPad = .012, displayGenomeBorder = F,
              genomeTicksLen = 3, genomeTicksTextSize = 0, genomeTicksScale = 5000000000,
              genomeLabelTextSize = 12, genomeLabelDy = -47, genomeTicksDisplay = F)
  })
}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)
