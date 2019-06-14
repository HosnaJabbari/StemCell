library(fcros)
library(plyr)
library(Biostrings)
library(EnhancedVolcano)
library(ComplexHeatmap)

#Haley Greenyer 
#Updated 6/4/2019

#import the file 
setwd("C:/Users/hgree/Documents/Senior Year/Research Project")
df <- read.csv("Hosna-copy2.csv", header = TRUE)

#isolate relavant columns (gene differentiations) from the data set 
col1 = df$X.HUV.iPS..normalized.
col2 = df$X.EC.Diff..normalized.
col3 = df$X.Nn.Diff..normalized.

#Calculate differences with reference to HUV.iPS (stem cell)
diff1 = col1 - col2 #HUV.iPS - EC.Diff (stem-endothelial)
diff2 = col1 - col3 #HUV.iPS - Nn.Diff (stem-neuronal)

#Assign these differences to variables in the dataframe 
df$Diff1 = diff1
df$Diff2 = diff2

#create subsets of the dataframe based on new variables 
#sub1 = "upregulated", sub2 = "downregulated" (in EC with respect to HUV)
sub1 <- subset(df, df$Diff1 < 0 & df$Diff2 > 0)
sub2 <- subset(df, df$Diff1 > 0 & df$Diff2 < 0)
sub3 = df

#add unique ID column to subsets for later identifcation 
ID1 <- 1:length(sub1$Probe.Set.ID); ID1
sub1$ID = ID1

ID2 <- 1:length(sub2$Probe.Set.ID); ID2
sub2$ID = ID2

#create new dataframes for fcros analysis (sub1.1 = "upregulated", sub2.1 = "downregulated" with respect to EC)
sub1.1 <- data.frame("HUV"= sub1$X.HUV.iPS..normalized., "EC" = sub1$X.EC.Diff..normalized., "Nn" = sub1$X.Nn.Diff..normalized., "HUVEC" = sub1$X.HUVEC..normalized., "HFN" = sub1$X.HFN..normalized., "ID" = sub1$Probe.Set.ID)

sub2.1 <- data.frame("HUV" = sub2$X.HUV.iPS..normalized., "EC" = sub2$X.EC.Diff..normalized.,"Nn" = sub2$X.Nn.Diff..normalized.,"HUVEC" = sub2$X.HUVEC..normalized., "HFN" = sub2$X.HFN..normalized., "ID" = sub2$Probe.Set.ID)

#define test condition labels (our two differentiated conditions and reference conditions)
test1 <- c("HUVEC","EC")
test2 <- c("HFN", "Nn")

#the induced stem cell condition ("HUV.iPS") is used as the control for all tests 
#(up and downregulation is with respect to HUV.iPS) 

fc1.1 <- fcros(sub1.1, "HUV", test1, 0, 1) #test for upregulation in EC
fc1.2 <- fcros(sub1.1, "HUV", test2, 0, 1) #test for downregulation in Nn

fc2.1 <- fcros(sub2.1, "HUV", test1, 0, 1) #test for downregulation in EC
fc2.2 <- fcros(sub2.1, "HUV", test2, 0, 1) #test for upregulation in Nn

#Add fcros characteristics to up and down-regulated subsets 
#fold changes 
sub1$FC.1 <- fc1.1$FC 
sub1$FC.2 <- fc1.2$FC
sub2$FC.1 <- fc2.1$FC
sub2$FC.2 <- fc2.2$FC

#f-value (adjusted p values)
sub1$f.val.1 <- fc1.1$f.value
sub1$f.val.2 <- fc1.2$f.value
sub2$f.val.1 <- fc2.1$f.value
sub2$f.val.2 <- fc2.2$f.value

#upregulation levels : FC > 2, f.val > 0.95 
#downregulation levels : FC < 0.5, f.val < 0.05 

#extract significant expression levels 
#subset 1 (upregulated in EC and downregulated in Nn)
sub1.EC.up <- subset(sub1, sub1$FC.1 >= abs(2) & sub1$f.val.1 >= 0.95) 
sub1.Nn.down <- subset(sub1, sub1$FC.2 <= abs(0.5) & sub1$f.val.2 <= 0.05)

#subset 2 (downregulated in EC and upregulated in Nn)
sub2.EC.down <- subset(sub2, sub2$FC.1 <= abs(0.5) & sub2$f.val.1 <= 0.05) 
sub2.Nn.up <- subset(sub2, sub2$FC.2 >= abs(2) & sub2$f.val.2 >= 0.95)

#extract ID's from each subest and determine intersection 
common1 <- intersect(sub1.EC.up$ID, sub1.Nn.down$ID)
common2 <- intersect(sub2.EC.down$ID, sub2.Nn.up$ID)

#select the intersection genes from sub1 and sub2 
sub1.sig <- sub1[common1,]
sub2.sig <- sub2[common2,]

#extract subsets for heatmaps 
#upregulated in EC and downregulated in Nn with respect to iPS
hsub1 <- as.matrix(sub1[ ,c(3:5)])
heatmap(hsub1,cexCol = 1)

#downregulated in EC and upregulated in Nn with respect to iPS
hsub2 <- as.matrix(sub2[ ,c(3:5)])
heatmap(hsub2, cexCol = 1)
