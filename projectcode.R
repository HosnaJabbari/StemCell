library(fcros)
library(plyr)
library(Biostrings)
library(EnhancedVolcano)
library(ComplexHeatmap)

#Haley Greenyer 
#Updated 6/19/2019

#import the file 
setwd("C:/Users/hgree/Documents/Senior Year/Research Project")
df <- read.csv("Hosna-copy2.csv", header = TRUE)

#isolate relavant columns (gene differentiations) from the data set 
col1 = df$X.HUV.iPS..normalized.
col2 = df$X.EC.Diff..normalized.
col3 = df$X.Nn.Diff..normalized.

#Calculate differences in expression levels with reference to HUV.iPS (stem cell)
diff1 = col1 - col2 #HUV.iPS - EC.Diff (stem-endothelial)
diff2 = col1 - col3 #HUV.iPS - Nn.Diff (stem-neuronal)

#Assign these differences to variables in the dataframe 
df$EC.Diff = diff1
df$Nn.Diff = diff2

#create subsets of the dataframe based on new variables 
#up regulated in EC and downregulated in Nn
sub1 <- subset(df, df$EC.Diff < 0 & df$Nn.Diff > 0)
#down regulated in EC and upregulated in Nn
sub2 <- subset(df, df$EC.Diff > 0 & df$Nn.Diff < 0)

#both upregulated
sub.up <- subset(df, df$EC.Diff < 0 & df$Nn.Diff < 0)
#both downregulated
sub.down <- subset(df, df$EC.Diff > 0 & df$Nn.Diff > 0)

#add unique ID column to subsets for later identifcation 
ID <- 1:length(df$Probe.Set.ID)
df$ID = ID

#the induced stem cell condition ("HUV.iPS") is used as the control for all fcros analyses
fcros.EC <- fcros(df, "X.HUV.iPS..normalized.", "X.EC.Diff..normalized.", 0, 1) #test for up regulation in EC
fcros.Nn <- fcros(df, "X.HUV.iPS..normalized.", "X.Nn.Diff..normalized.", 0, 1) #test for down regulation in Nn

#Add fcros characteristics to the data frame 
df$FC.EC <- fcros.EC$FC #fold changes 
df$FC.Nn <- fcros.Nn$FC
df$fval.EC <- fcros.EC$f.value #f value (adjusted p value)
df$fval.Nn <- fcros.Nn$f.value

#log2 scale for fold changes 
#upregulation levels : FC > 2, f.val > 0.95 
#downregulation levels : FC < 0.5, f.val < 0.05 

#extract significantly downregulated genes from fold changes 
EC.sig.up <- subset(df, df$FC.EC >= 2) 
Nn.sig.up <- subset(df, df$FC.Nn >= 2)

#extract significanatly downregulated genes from fold changes 
EC.sig.down <- subset(df, df$FC.EC <= 0.5)
Nn.sig.down <- subset(df, df$FC.Nn <= 0.5)

#create subsets that include the f value significance levels
EC.up.fvalue <- subset(df, df$FC.EC >= 2 & df$fval.EC >=0.95) 
Nn.up.fvalue <- subset(df, df$FC.Nn >= 2 & df$fval.Nn >=0.95)
EC.dwn.fvalue <- subset(df, df$FC.EC <= 0.5 & df$fval.EC <= 0.05)
Nn.dwn.fvalue <- subset(df, df$FC.Nn <= 0.5 & df$fval.Nn <=0.05)

#Determine desired intersections 
common1 <- intersect(EC.sig.up$ID, Nn.sig.down$ID) #up regulated in EC and down regulated in Nn
common2 <- intersect(EC.sig.down$ID, Nn.sig.up$ID) #down regulated in EC and up regulated in Nn
common3 <- intersect(EC.sig.up$ID, Nn.sig.up$ID) #up regulated in both 
common4 <- intersect(EC.sig.down$ID, Nn.sig.down$ID) #down regulated in both 

#determine intersections considering f value 
common5 <- intersect(EC.up.fvalue$ID, Nn.dwn.fvalue$ID)
common6 <- intersect(EC.dwn.fvalue$ID, Nn.up.fvalue$ID)
common7 <- intersect(EC.up.fvalue$ID, Nn.up.fvalue$ID)
common8 <- intersect(EC.dwn.fvalue$ID, Nn.dwn.fvalue$ID)

#select the intersection genes from the data frame to make desired subsets 
EC.up.Nn.dwn <- df[common1,]
EC.dwn.Nn.up <- df[common2,]
EC.Nn.up <- df[common3,]
EC.Nn.dwn <- df[common4,]

#select the intersection genes from the data frame considering f value 
EC.up.Nn.dwn.fval <- df[common5,]
EC.dwn.Nn.up.fval <- df[common6,]
EC.Nn.up.fval <- df[common7,]
EC.Nn.dwn.fval <- df[common8,]

#extract subsets for heatmaps 
#up regulated in EC and down regulated in Nn with respect to HUV.iPS
hsub1 <- as.matrix(sub1[ ,c(3:5)])
heatmap(hsub1,cexCol = 1)

#downregulated in EC and upregulated in Nn with respect to HUV.iPS
hsub2 <- as.matrix(sub2[ ,c(3:5)])
heatmap(hsub2, cexCol = 1)

#up regulatwed in both with respect to HUV.iPS
hsub3 <- as.matrix(sub.up[ ,c(3:5)])
heatmap(hsub3,cexCol = 1)

#down regulated in both with resoect to HUV.iPS
hsub4 <- as.matrix(sub.up[ ,c(3:5)])
heatmap(hsub4,cexCol = 1)

#Sort Gene Ontology subsets for each gene expression level set 
#Transcription factors and cofactors 
trans.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("transcription", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
trans.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("transcription", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
trans.up <- EC.Nn.up[grep("transcrption", EC.Nn.up$Gene.Ontology.Biological.Process),]
trans.dwn <- EC.Nn.dwn[grep("transcription", EC.Nn.dwn$Gene.Ontology.Biological.Process),]

#epigenetic regulators 
epig.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("chromatin", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
epig.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("chromatin", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
epig.up <- EC.Nn.up[grep("chromatin", EC.Nn.up$Gene.Ontology.Biological.Process),]
epig.dwn <- EC.Nn.dwn[grep("chromatin", EC.Nn.dwn$Gene.Ontology.Biological.Process),]

#signaling 
signal.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("signal", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
signal.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("signal", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
signal.up <- EC.Nn.up[grep("signal", EC.Nn.up$Gene.Ontology.Biological.Process),]
signal.dwn <- EC.Nn.dwn[grep("signal", EC.Nn.dwn$Gene.Ontology.Biological.Process),]

#metabolism 
met.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("metab", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
met.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("metab", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
met.up <- EC.Nn.up[grep("metab", EC.Nn.up$Gene.Ontology.Biological.Process),]
met.dwn <- EC.Nn.dwn[grep("metab", EC.Nn.dwn$Gene.Ontology.Biological.Process),]

#cell adhesison 
adh.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("adhesion", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
adh.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("adhesion", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
adh.up <- EC.Nn.up[grep("adhesion", EC.Nn.up$Gene.Ontology.Biological.Process),]
adh.dwn <- EC.Nn.dwn[grep("adhesion", EC.Nn.dwn$Gene.Ontology.Biological.Process),]

#extracellular matrix protein 
ext.EC.up.Nn.dwn <- EC.up.Nn.dwn[grep("extracellular", EC.up.Nn.dwn$Gene.Ontology.Biological.Process),]
ext.EC.dwn.Nn.up <- EC.dwn.Nn.up[grep("extracellular", EC.dwn.Nn.up$Gene.Ontology.Biological.Process),]
ext.up <- EC.Nn.up[grep("extracellular", EC.Nn.up$Gene.Ontology.Biological.Process),]
ext.dwn <- EC.Nn.dwn[grep("extracellular", EC.Nn.dwn$Gene.Ontology.Biological.Process),]
