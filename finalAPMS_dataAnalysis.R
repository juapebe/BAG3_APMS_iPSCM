####LOAD DEPENDENCIES####
library(reshape2)
library(ggplot2)
library(MSstats)
library(pheatmap)
source("./process4MSStats_3.R")
source("./auxiliaryFunctions.R")

####LOAD AND PREPARE FILES####
#Annotation is a list with of conditions and replicates
annotation<-read.table("./SetA+B+C+D_Fusion/annotation.txt", sep="\t", 
                  header=T, stringsAsFactors = F) #file is made from the make_annotation custom function.

#Uniprot database. for uniprotID<>geneName
uniprot<-read.delim("./referenceFiles/uniprot_singleID.txt",
                    header = TRUE, stringsAsFactors = F)[-1]

#An in-house list of genes associated to DCM. Compiled by Luke Judge
dcmSetAll<-read.delim("./referenceFiles/CardiomyopathyGenes_ExpandedPanel.txt")
dcmSet<-dcmSetAll[c("Combined", "UniprotID")] #This is the set that was used for the venn diagram

#Custom script for parsing maxquant files
l <- parse_maxquant(evidenceFile = "./SetA+B+C+D_Fusion/txt/evidence.txt", 
               proteinGroupsFile = "./SetA+B+C+D_Fusion/txt/proteinGroups.txt", 
               annotationFile = "./SetA+B+C+D_Fusion/annotation.txt")
#Data frame of intensities for PEPTIDES per EXPERIMENT. 'wide' format
peptideWide <- l[[2]]
peptideLong <- melt(data = peptideWide, value.var="Intensity")
peptideLong <- merge(x = peptideLong, y = annotation, by.x="variable", by.y="BioReplicaSaint", all.x = T)
peptideLong <- peptideLong[,-which(names(peptideLong) %in% c("Raw.file", "IsotopeLabelType", "SAINT"))]
peptideLong <- data.frame(run = peptideLong$variable, protein = peptideLong$Proteins, 
                          peptideSequenceCharge=peptideLong$feature, condition=peptideLong$Condition,
                          replicate = peptideLong$BioReplicate, runOrder = peptideLong$Run,
                          intensity=peptideLong$value)
peptideLong <- peptideLong[order(peptideLong$intensity, decreasing=T),]
peptideLong <- peptideLong[order(peptideLong$runOrder, decreasing=F),]

#Intensities for PROTEINS per EXPERIMENT. 'wide' format
proteinWide <- l[[1]]
proteinLong <- melt(data = proteinWide, value.var="Intensity")
proteinLong <- merge(x = proteinLong, y = annotation, by.x="variable", by.y="BioReplicaSaint", all.x = T)
proteinLong <- proteinLong[,-which(names(proteinLong) %in% c("Raw.file", "IsotopeLabelType", "SAINT"))]
proteinLong <- data.frame(run = proteinLong$variable, protein = proteinLong$Proteins, 
                          condition=proteinLong$Condition,
                          replicate = proteinLong$BioReplicate, runOrder = proteinLong$Run,
                          intensity=proteinLong$value)
proteinLong <- proteinLong[order(proteinLong$intensity, decreasing=T),]
proteinLong <- proteinLong[order(proteinLong$runOrder, decreasing=F),]

#List of removed contaminants 
contaminants <- l[[3]]

#Remove all samples we are not interested in
proteinLong <- proteinLong[proteinLong$condition %in% c("WTc11", "A4", "K22"),]
proteinLong$condition <- droplevels(proteinLong$condition)
peptideLong <- peptideLong[peptideLong$condition %in% c("WTc11", "A4", "K22"),]
peptideLong$condition <- droplevels(peptideLong$condition)

###QC plots
#Intensity of the iRT peptides per run
#most common charge is 2. use only that one for normalization and plotting
irtPeptides <- peptideLong[peptideLong$protein=="Biognosys|iRT-Kit_WR_fusion" & grepl("__2", peptideLong$peptideSequenceCharge),]
irtPeptides <- irtPeptides[!grepl("(ac)", irtPeptides$peptideSequenceCharge),]#remove acetylated peptides
p <- ggplot(data = irtPeptides, 
            mapping = aes(x=run, y=log2(intensity), color=peptideSequenceCharge, group = peptideSequenceCharge))
p+geom_point()+ ylim(c(0, 32)) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = "./plots/iRT_peptide_intensities.pdf", device = "pdf")

#Intensity of the BAG3 peptides
bag3Peptides <- peptideLong[peptideLong$protein=="O95817",]
p <- ggplot(data = bag3Peptides, 
            mapping = aes(x=as.character(sprintf("%02i",runOrder)), y=log2(intensity), color=peptideSequenceCharge, group = peptideSequenceCharge))
p+geom_point() + geom_line() +  ylim(c(0, 32)) + theme(legend.position="none") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = "./plots/bag3Peptide_intensitites.pdf", device = "pdf")

