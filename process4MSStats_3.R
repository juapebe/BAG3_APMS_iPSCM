#REMEMBER TO DEAL WITH BAG3 VARIANTS!

library(MSstats)
library(ggplot2)
library(reshape2)
library(pheatmap)


####UTILIZING SAINTQ OUTPUTS####
#Load Saintq output files
saintqScores<-read.delim("./SetA+B+C+D_Fusion/saintq/", header=T, stringsAsFactors = F)
saintqScores[grepl("O95817",saintqScores$Prey),'Prey']<-'O95817'#To deal with BAG3 variants

cutoff<-0.1 #Cutoff for the FDR

#When passing to matrix form, using fun.aggregate=min deals with multiple values for BAG3,
#while keeping missing values as 1 (for fillers)
bfdrMatrix<-dcast(saintqScores, formula = Prey ~ Bait, value.var = "BFDR", drop=F, fun.aggregate=min, na.rm=T, fill = 0) #casting to matrix form
peptMatrix<-dcast(saintqScores, formula = Prey ~ Bait, value.var = "X.Pep", drop=F, fun.aggregate=min, na.rm=T, fill = 0) #casting to matrix form
row.names(bfdrMatrix)<-bfdrMatrix[,1]; bfdrMatrix<-bfdrMatrix[,-which(names(bfdrMatrix) %in% c("Prey"))]
row.names(peptMatrix)<-peptMatrix[,1]; peptMatrix<-peptMatrix[,-which(names(peptMatrix) %in% c("Prey"))]

#How many significant interactors for each prey?
saintqScores$Significant<-saintqScores$BFDR<=cutoff
png(paste("./#SignificantHits_FDR", cutoff, ".png", sep=""))
p<-ggplot(saintqScores, aes(x=Bait, fill=factor(Significant)))
p+geom_bar()+labs(title=paste("# Significant hits with FDR below", cutoff), x="Bait") #plotting it
dev.off()
sigTable<-apply(bfdrMatrix, MARGIN = 2, function(x) table(x<=cutoff)) #In table form
write.table(sigTable, paste("./#SignificantHits_FDR", cutoff, ".txt", sep=""), quote=F, sep="\t")


#A matrix with all interactors that are significant for at least one of the baits
oneSignificant<-apply(bfdrMatrix, MARGIN = 1, function(x) any(x<=0.05, na.rm = T))#the names
bfdrMatrixSig<-bfdrMatrix[oneSignificant==TRUE,]

#a subset that contains entries depending on only 1 or 2 peptides (to be handled carefully)
#note: I do those that contain 1 or 2 entries in all of them - otherwise they are kept normally
peptMatrixSig<-peptMatrix[row.names(peptMatrix) %in% row.names(bfdrMatrixSig),]
lowPept<-apply(peptMatrixSig, MARGIN = 1, function(x) all(x<=2, na.rm = T))#the names
bfdrLowPeptMatrix<-bfdrMatrixSig[lowPept==TRUE,]
lowPeptMatrixSig<-peptMatrixSig[lowPept==TRUE,]




##BFDRLOWPEPTMATRIX IS THE ONE WITH SIGNIFICANT PROTEINS IDENTIFIED BY ONLY FEW PEPTIDES
##BFDRMATRIXSIG IS THE ONE WITH SIGNIFICANT PROTEINS (ALL)


####MSSTATS PROCESSING####

#Setup Files
#NOTE: using the ones from the 
evidence<-read.table("./SetA+B+C+D_Fusion/txt/evidence.txt", sep="\t", 
                    header=T, stringsAsFactors = F)

annot<-read.table("./SetA+B+C+D_Fusion/annotation.txt", sep="\t", 
                  header=T, stringsAsFactors = F) #file is made from the make_annotation custom function.

#Using "annot2" every sample is treated separately.
#annot2 <- annot
#annot2$Condition <- paste(annot2$Condition, "_#", 1:nrow(annot2), sep="")

proteinGroups<-read.table("./SetA+B+C+D_Fusion/txt/proteinGroups.txt", sep="\t", 
                           header=T, stringsAsFactors = F)

#Uniprot database. for naming et al
uniprot<-read.delim("~/Protein_Databases/Uniprot_human_Jul2016/uniprot_singleID.txt",
                    header = TRUE, stringsAsFactors = F)


#Filtering (or not) the evidence entries to those associated with SAINTq significant hits
# infileSig<-infile2[infile2$Leading.razor.protein %in% row.names(bfdrMatrixSig),]
# print(paste("Number of unique proteins in total: ", length(unique(infile2$Leading.razor.protein))))
# print(paste("Number of unique *significant* proteins in total: ", length(unique(infileSig$Leading.razor.protein))))


msStatsInput<-MaxQtoMSstatsFormat(evidence = evidence, annotation = annot, 
                                  proteinGroups = proteinGroups)

#OPTIONAL: narrow down to only CM samples of interest
msStatsInput <- msStatsInput[!grepl("iPS_|HSPB7|RBM20", msStatsInput$Condition),]

# msStatsInputSig<-MaxQtoMSstatsFormat(evidence = infileSig, annotation = annot, 
#                                      proteinGroups = proteinGroups)


####MSSTATS PROCESSING####

###STEP 1 - dataProcess
#Run DataProcess. All hits
#msStatsProc_medianNorm<-dataProcess(msStatsInput) #default options
msStatsProc_noNorm<-dataProcess(msStatsInput, normalization = FALSE,
                                cutoffCensored = "minFeature", summaryMethod = "TMP", 
                                equalFeatureVar = FALSE, MBimpute = TRUE,
                                fillIncompleteRows = TRUE)

msStatsProc_medianNorm<-dataProcess(msStatsInput, normalization = 'equalizeMedians', 
                                  cutoffCensored = "minFeature", summaryMethod = "TMP",
                                  equalFeatureVar = TRUE, MBimpute = TRUE,
                                  fillIncompleteRows = TRUE)

#For normalization on global standards (BAG3), limit it to the ions that are (a) not present in the control (one rep has a few there) 
#and (b) present in all others (NOTE: this excluded by now) and (c) not phosphorylated (NOTE: not used, as it seemed to throw off normalization)
bag3Peptides <- msStatsInput[msStatsInput$ProteinName=="O95817" & !is.na(msStatsInput$Intensity),]
bag3PeptidesPhosphorylated <- bag3Peptides[grepl(pattern = "\\(ph\\)", x = bag3Peptides$PeptideSequence),]
bag3PeptidesPhosphorylated <- gsub(pattern = "\\(ph\\)", replacement = "", x = bag3PeptidesPhosphorylated$PeptideSequence)#only the sequence
bag3PeptidesInControl <- bag3Peptides[bag3Peptides$Condition=="WTc11","PeptideSequence"]
bag3Peptides <- bag3Peptides[!bag3Peptides$PeptideSequence %in% bag3PeptidesInControl,]#removes those in control
bag3Peptides <- bag3Peptides[!bag3Peptides$PeptideSequence %in% bag3PeptidesPhosphorylated,]#removes those with phosphorylation
bag3Peptides$PeptideSequenceCharge <- paste(bag3Peptides$PeptideSequence, bag3Peptides$PrecursorCharge, sep="_")
bag3PeptideStandards <- as.vector(unique(droplevels(bag3Peptides$PeptideSequence)))

msStatsProc_baitNorm<-dataProcess(msStatsInput, normalization = 'globalStandards', nameStandards = bag3PeptideStandards, 
                                    cutoffCensored = "minFeature", summaryMethod = "TMP",
                                    equalFeatureVar = TRUE, MBimpute = TRUE,
                                    fillIncompleteRows = TRUE)

#Run DataProcess. Only significant hits
# msStatsProcSig_noNorm<-dataProcess(msStatsInputSig, normalization = FALSE,
#                                    nameStandards="O95817", cutoffCensored = "minFeature",
#                                    summaryMethod = "TMP", equalFeatureVar = FALSE, MBimpute = TRUE,
#                                    fillIncompleteRows = TRUE)
# msStatsProcSig_baitNorm<-dataProcess(msStatsInputSig, normalization = 'globalStandards', 
#                                      nameStandards="O95817", cutoffCensored = "minFeature",
#                                      summaryMethod = "TMP", equalFeatureVar = FALSE, MBimpute = TRUE,
#                                      fillIncompleteRows = TRUE)


#Plots - all hits
dir.create("./msStats/", showWarnings = TRUE)
dataProcessPlots(msStatsProc_noNorm, type = "ProfilePlot", address = "./msStats/BAG3_noNormalized", 
                 which.Protein = "O95817", save_condition_plot_result = T, text.size=2) #per protein Plot with intensities
dataProcessPlots(msStatsProc_noNorm, type = "QCPlot", address = "./msStats/BAG3_noNormalized", 
                 text.size = 2, which.Protein = "O95817") #Quality control Plot

dir.create("./msStats/", showWarnings = TRUE)
dataProcessPlots(msStatsProc_medianNorm, type = "ProfilePlot", address = "./msStats/BAG3_medianNorm", 
                 which.Protein = "O95817", save_condition_plot_result = T, text.size=2) #per protein Plot with intensities
dataProcessPlots(msStatsProc_medianNorm, type = "QCPlot", address = "./msStats/BAG3_medianNorm", 
                 which.Protein = "O95817", text.size = 2) #Quality control Plot

dataProcessPlots(msStatsProc_baitNorm, type = "ProfilePlot", address = "./msStats/BAG3_baitNorm", 
                 which.Protein = "O95817", save_condition_plot_result = T, text.size=2) #per protein Plot with intensities
dataProcessPlots(msStatsProc_baitNorm, type = "QCPlot", address = "./msStats/BAG3_baitNorm", 
                 which.Protein = "O95817", text.size = 2) #Quality control Plot

# #Plots - only significant hits
# dir.create("./plots_sig_noNorm", showWarnings = TRUE)
# dataProcessPlots(msStatsProcSig_noNorm, type = "ProfilePlot", address = "./plots_sig_noNorm/") #per protein Plot with intensities
# dataProcessPlots(msStatsProcSig_noNorm, type = "QCPlot", address = "./plots_sig_noNorm/") #Quality control Plot
# 
# dir.create("./plots_sig_baitNorm", showWarnings = TRUE)
# dataProcessPlots(msStatsProcSig_baitNorm, type = "ProfilePlot", address = "./plots_sig_baitNorm/") #per protein Plot with intensities
# dataProcessPlots(msStatsProcSig_baitNorm, type = "QCPlot", address = "./plots_sig_baitNorm/") #Quality control Plot



###STEP 1.5 - QC of samples
#Clustering of samples (NOTE that this is only for the 'significant' hits from SAINTq!)
sampledf<-msStatsProc_medianNorm$ProcessedData
sampledf$SAMPLE<-paste(sampledf$GROUP_ORIGINAL, sampledf$RUN, sep = "_")
sampleMatrix<-dcast(data = sampledf, FEATURE ~ SAMPLE, 
                    value.var = "ABUNDANCE", drop = FALSE)
rownames(sampleMatrix)<-sampleMatrix$FEATURE
sampleMatrix<-sampleMatrix[-1]
sampleMatrix_cor <- cor(sampleMatrix, use="pairwise.complete.obs", method="pearson")
pheatmap(sampleMatrix_cor)


#Checking BAG3 total counts. NOTE: QC plots from MSSTATS are better!
sampledfBAG3<-sampledf[sampledf$PROTEIN=="O95817",]
sampledfBAG3$PROTEIN<-factor(sampledfBAG3$PROTEIN)
sampledfBAG3$FEATURE<-factor(sampledfBAG3$FEATURE)
sampleMatrixBAG3<-dcast(data = sampledfBAG3, PROTEIN ~ SAMPLE, 
                        value.var = "ABUNDANCE", drop = FALSE, fun.aggregate = sum, na.rm=T)
rownames(sampleMatrixBAG3)<-sampleMatrixBAG3$PROTEIN
sampleMatrixBAG3<-t(sampleMatrixBAG3[-1])
sampleMatrixBAG3<-as.data.frame(sampleMatrixBAG3)
#sampleMatrixBAG3_cor <- cor(sampleMatrixBAG3, use="pairwise.complete.obs", method="pearson")
g<-ggplot(data = sampleMatrixBAG3, aes(x=row.names(sampleMatrixBAG3), y=O95817))+
    geom_bar(stat="Identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
g


###STEP 2 - GROUP COMPARISONS
#checkLevels
levels(msStatsProc_noNorm$ProcessedData$GROUP_ORIGINAL)

#Define which comparisons I want
comparisonTemplate <- rep(0, length(levels(msStatsProc_noNorm$ProcessedData$GROUP_ORIGINAL)))
comparison_A4vsCMcontrol <- matrix(c(1, 0, 0, 0, 0, 0, 0, -1), nrow=1)
comparison_homoDMSOvsCMcontrol <- matrix(c(0, 0, 0, 1, 0, 0, 0, -1), nrow=1)
comparison_homoDMSOvsTetOnBAG3 <- matrix(c(0, 0, 0, 1, 0, 0, -1, 0), nrow=1)

comparison<-rbind(comparison_A4vsCMcontrol, comparison_homoDMSOvsCMcontrol, comparison_homoDMSOvsTetOnBAG3)
row.names(comparison)<- c("comparison_A4vsCMcontrol", "comparison_homoDMSOvsCMcontrol", "comparison_homoDMSOvsTetOnBAG3")
msStatsComparison<-groupComparison(comparison, data=msStatsProc_noNorm)


###STEP 2.5 - QC of comparisons
#Volcano Plot
#compDataFilt<-compData[is.na(compData$issue),]

groupComparisonPlots(data= msStatsComparison$ComparisonResult, type="VolcanoPlot", address = "VolcanoPlotTest.pdf", sig = 0.2)

#Heat Map
groupComparisonPlots(data=msStatsComparison$ComparisonResult,  type="Heatmap")
