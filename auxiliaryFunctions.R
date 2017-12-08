make_annotation<- function(expDesignFiles, controlNames, outputName){
#A function to make the (keys) file that is used by both SAINTq('keys') and MsStats
#NOTE: good practice to have your ExpDesignTemplate to include a 'Condition' column.
  #MODIFIED: now concatenates multiple annotation files
  #controlNames is a vector with the string names of the control conditions
  #e.g. c("GFP", "Vector")
  expDesignList <- lapply(expDesignFiles, read.table, sep="\t", header=T, stringsAsFactors = F)
  print(expDesignList)
  expDesign <- do.call(rbind, expDesignList)
  expDesign=expDesign[expDesign$Experiment!="bsa" & expDesign$Experiment!="yeast",] #remove blanks
  
  BioReplicaSaint<-matrix(ncol = 2, nrow = dim(expDesign[2]))
  for(c in 1:length(expDesign$Condition)){
    BioReplicaSaint[c,] <- c(expDesign$Condition[c], sum(BioReplicaSaint[,1]==expDesign$Condition[c], na.rm = T)+1)
  }
  BioReplicaSaint2 = paste(BioReplicaSaint[,1], BioReplicaSaint[,2], sep="_")
  print(head(expDesign))
  f<-data.frame(Raw.file=expDesign$Name, Condition = expDesign$Condition,
                BioReplicate = BioReplicaSaint[,2], Run=1:length(expDesign$Name), 
                IsotopeLabelType="L", SAINT="T", BioReplicaSaint=BioReplicaSaint2, 
                stringsAsFactors = F)
  #sets the controls
  f[f$Condition %in% controlNames,]$SAINT<-"C"
  write.table(f, file = outputName, sep="\t", col.names=T, row.names=F, quote = F)
}

removeLowPept <- function(m, proteinGroups=proteinGroups){
#adds column to peptide matrix, marking those from low peptide proteins
  proteinGroupsShort <- proteinGroups
  proteinGroupsShort$Protein.IDs <- gsub(";.*", "", proteinGroupsShort$Protein.IDs)
  proteinGroupsShort$Peptide.counts..razor.unique. <- as.numeric(gsub(";.*", "", proteinGroupsShort$Peptide.counts..razor.unique.))
  print(str(proteinGroupsShort))
  proteinGroupsShort <- proteinGroupsShort[proteinGroupsShort$Peptide.counts..razor.unique.<2,]
  
  m[row.names(m) %in% proteinGroupsShort$Protein.IDs,1] <- "LOWPEPT"
  return(m)
}

impute_missing_values <- function(v){
  #Using the "Perseus method": inpute missing as a random sample from a normal distribution 2.5 stdevs
  #below actual nonzero population, and an SD of 0.3*population
  plot(density(log10(v[v!=0])))
  m <- mean(log10(v[v!=0]))
  sd <- sd(log10(v[v!=0]))
  v[v==0] <- 10^(rnorm(n = length(v[v==0]), mean = m-2.5*sd, sd=0.3*sd)) #Using perseus defaults
  lines(density(log10(v)))
  return(v)
}

computeZ <- function(column){
#Computing Z-scores (only when taking the Inf values out of the equation)
  avg <- mean(column)
  stdev <- sd(column)
  #z <- (column-avg)/stdev
  pval <- pnorm(column, mean = avg, sd = stdev, lower.tail = FALSE)
  return(pval)
}