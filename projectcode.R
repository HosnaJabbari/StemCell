library(fcros)
library(plyr)
library(Biostrings)
library(EnhancedVolcano)
library(ComplexHeatmap)

#Haley Greenyer 
#Updated 5/18/2019

#import the file 
setwd("C:/Users/hgree/Documents/Senior Year/Research Project")
df <- read.csv("Hosna-copy2.csv", header = TRUE)

#isolate relavant columns from the data set 
col1 = df$X.HUV.iPS..normalized.
col2 = df$X.EC.Diff..normalized.
col3 = df$X.Nn.Diff..normalized.

#Calculate differences with reference to HUV.iPS
diff1 = col1 - col2 #HUV.iPS - EC.Diff
diff2 = col1 - col3 #HUV.iPS - Nn.Diff 

#Assign these differences to variables in the dataframe 
df$Diff1 = diff1
df$Diff2 = diff2

#create subsets of the dataframe based on new variables 
#sub1 = "upregulated", sub2 = "downregulated" (in EC)
sub1 <- subset(df, df$Diff1 < 0 & df$Diff2 > 0)
sub2 <- subset(df, df$Diff1 > 0 & df$Diff2 < 0)

#add unique ID column to subsets for later identifcation 
ID1 <- 1:length(sub1$Probe.Set.ID); ID1
sub1$ID = ID1

ID2 <- 1:length(sub2$Probe.Set.ID); ID2
sub2$ID = ID2

#create new dataframes for fcros analysis (sub1.1 = "upregulated", sub2.1 = "downregulated" with respect to EC)
sub1.1 <- data.frame("HUV"= sub1$X.HUV.iPS..normalized., "EC" = sub1$X.EC.Diff..normalized., "Nn" = sub1$X.Nn.Diff..normalized., "ID" = sub1$Probe.Set.ID)

sub2.1 <- data.frame("HUV" = sub2$X.HUV.iPS..normalized., "EC" = sub2$X.EC.Diff..normalized.,"Nn" = sub2$X.Nn.Diff..normalized., "ID" = sub2$Probe.Set.ID)

#define test condition labels 
test <- c("EC","Nn")

#perform fcros for subsets, where fc1 is upregulated in EC and downregulated in Nn, fc2 is opposite 
fc1 <- fcros(sub1.1, "HUV", test, 0, 1)
fc2 <- fcros(sub2.1, "HUV", test, 0, 1)

#fcros for entire set
fsub <- data.frame("HUV" = df$X.HUV.iPS..normalized.,"EC" = df$X.EC.Diff..normalized.,"Nn" = df$X.Nn.Diff..normalized.)
fc3 <-fcros(fsub, "HUV", test, 0, 1)

#select top 1000 down/up regulated genes of all subsets 
T.fc1 <- fcrosTopN(fc1,1000)
T.fc2 <- fcrosTopN(fc2,1000)
T.fc3 <- fcrosTopN(fc3,1000)

#Add fcros characteristics to up and down-regulated subsets 
sub1$FC <- fc1$FC
sub2$FC <- fc2$FC
sub1$p.val <- fc1$p.value
sub2$p.val <- fc2$p.value

#extract significantly significant subsets to a list (5% significance level)
sub1.sig <- subset(sub1, sub1$FC > abs(2) & sub1$p.val< 0.05)
sub2.sig <- subset(sub2, sub2$FC > abs(2) & sub2$p.val< 0.05)
 
#volcano plots w/ optional labels 
#upregulated (in EC)
plot(fc1$FC, -log10(fc1$p.value), pch=20, main="Volcano plot")
plot(fc1$FC, -log10(fc1$p.value), pch=20, col=ifelse(fc1$p.value < 0.05 & fc1$FC > abs(2), "red", "black"), xlim = c(0,20))
#text(sub1.sig$FC, -log10(sub1.sig$p.val),labels = sub1.sig$Gene.Symbol)

#downregulated (in EC)
plot(fc2$FC, -log10(fc2$p.value), pch=20, main="Volcano plot")
plot(fc2$FC, -log10(fc2$p.value), pch=20, col=ifelse(fc2$p.value < 0.05 & fc2$FC > abs(2), "red", "black"))
#text(sub2.sig$FC, -log10(sub2.sig$p.val),labels = sub2.sig$Gene.Symbol)

#extract subsets for heatmaps 
#upregulated in EC and downregulated in Nn with respect to iPS
hsub1 <- as.matrix(sub1[ ,c(3:5)])
heatmap(hsub1)

#downregulated in EC and upregulated in Nn with respect to iPS
hsub2 <- as.matrix(sub2[ ,c(3:5)])
heatmap(hsub2)
