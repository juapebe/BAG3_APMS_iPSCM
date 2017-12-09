#! /usr/bin/Rscript --vanilla
#
# --------------------
# process4saintq.R
# --------------------
# 
# Description
# 
# Take the evidence file and keys (with the BioReplicaSaint column) 
# and it generates the input for saintq.
# --------------------
# 
# For any issue, contact David Jimenez-Morales at Krogan Lab
# 
# (C) developers @ Kroganlab
# 


#NOTE: make_annotation() file is now in the auxiliaryFunctions.R script.

######
#MOST OF THIS IS BASED ON DAVID JIMENEZ-MORALES' SCRIPT, PROCESS4SAINTQ.R. that is the original.
#Modified to be contained within a function that returns global (<<-) variables. 
#This is done before producing the SAINTq-processed objects.
parse_maxquant <- function(evidenceFile, proteinGroupsFile, annotationFile){
    #Load dependencies
    suppressMessages(require(bit64))
    suppressMessages(library(reshape2))
    suppressMessages(library(data.table))

    #NOTE: ran everything together for simplicity.
    data <- read.delim(evidenceFile, sep='\t', stringsAsFactors = F) #evidence.txt
    proteinGroups <- read.delim(proteinGroupsFile, sep="\t", stringsAsFactors = F) #proteinGroups.txt
    keys <- read.delim(annotationFile, sep='\t', stringsAsFactors = F) #from make_annotation
    
    ####
    #INSERTED BY JUAN rename Raw.file if instead of RawFile (consequence of different versions of MaxQuant)
    if(("Raw.file" %in% names(keys)) & !("RawFile" %in% names(keys))) names(keys)[names(keys)=="Raw.file"]<-"RawFile"
    if(("Raw.file" %in% names(data)) & !("RawFile" %in% names(data))) names(data)[names(data)=="Raw.file"]<-"RawFile"
    ####

    # Check that the keys file is correct
    if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run','SAINT','BioReplicaSaint') %in% colnames(keys))){
        cat('\nCOLNAMES IN KEYS NOT CONFORM TO SCHEMA. One of these is lost\n
            \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\\ntRun\n\n\tSAINT\tBioReplicaSaint\n\n')
        stop('Try again once revised\n')
    }

    # MERGING THE DATA
    # Checking that the keys make sense
    unique_data <- unique(data$RawFile)
    unique_keys <- unique(keys$RawFile)
    
    # Rawfiles on Keys not found on the data
    keys_not_found <- setdiff(unique_keys, unique_data)
    # Rawfiles on Data not found on the keys
    data_not_found <- setdiff(unique_data, unique_keys)
    
    if ( (length(keys_not_found) != 0) & ( length(data_not_found) != 0) ) {
        cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
        cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
        stop('\nThis script is sorry, but it needs to stop this because something is going on between your keys and evidence files so you better check\n')
    }
    
    ## select only required attributes from MQ format
    datamerged <- merge(data, keys, by='RawFile')
    
    #REPLACED BY JUAN FOR THIS
    data_f1 <- datamerged[grep(";",datamerged$Protein.group.IDs, invert=T),]
    ####

    # Get rid of the contaminants
    contaminants <- unique(data_f1[grep("CON__|REV__", data_f1$Proteins),c("Proteins", "Gene.names")])
    print("REMEMBER to check the object 'contaminants' in case there's something you didn't want to remove!")
    data_f2 <- data_f1[grep("CON__|REV__", data_f1$Proteins, invert=T),]
    
    if(length(which(data_f2$Proteins==""))>0)  data_f2 <- data_f2[-which(data_f2$Proteins==""),]
    
    
    # Set the intensity as numeric to avoid overflow problems
    data_f2$Intensity = as.numeric(data_f2$Intensity)


    # ADDED BY JUAN - remove Methionine-containing peptides (NOTE: should be made optional argument?)
    #data_f2<-data_f2[grep("M", data_f2$Sequence,  invert=T),]


    #MODIFIED BY JUAN - using protein groups as unique identifiers
    suppressWarnings(protBiorepIntensity  <- dcast(data=data_f2[,c("Leading.razor.protein","BioReplicaSaint","Intensity")], 
                                                   Leading.razor.protein~BioReplicaSaint, value.var = "Intensity", fun.aggregate = sum,
                                                   na.rm = T, fill=0))

    print("NOTE that the peptides are associated with the first protein in the protein group (for the sake of having single identifiers)")
    #Choose 'max' intensity instead of 'sum' to deal with duplicate entries (as MSstats does) 
    #(NOTE: implement as option? Im not sure which one works best) NEEDS TO BE TESTED
    data_f2<-data_f2[!is.na(data_f2$Intensity),] #remove NA Intensities
    data_f2$feature<-paste(data_f2$Modified.sequence, data_f2$Charge, sep = "_")
    peptideBiorepIntensity <- dcast(data=data_f2[,c("Leading.razor.protein","feature","BioReplicaSaint","Intensity")],
                                    Leading.razor.protein+feature~BioReplicaSaint, value.var = "Intensity", fun.aggregate = max, 
                                    na.rm = T, fill=0 )

    #distribution of frequencies of feature+replicate combinations (should be mostly '1')
    #table(table(paste(data_f2$feature, data_f2$BioReplicaSaint, sep = "/")))

    names(protBiorepIntensity)[1] <- "Proteins"; names(peptideBiorepIntensity)[1] <- "Proteins"
    protBiorepIntensity <- protBiorepIntensity[complete.cases(protBiorepIntensity),]
    # Just in case there is NA values:
    if(nrow(protBiorepIntensity[!complete.cases(protBiorepIntensity),])!=0){
        print("WARNING: there were some NA values in the list of intensities. is this normal?")
    }

    # Take only unique values of both sequences and proteins. Aggregate them to get the sequences for each protein
    data_f2uniques <- unique(data_f2[,c("Proteins", "Sequence")])
    protPep <- aggregate(Sequence ~ Proteins, data_f2uniques, FUN = paste, collapse="|" )
    

    # PROTEINS: extra step to add the information about peptides
    almost <- merge(protBiorepIntensity, protPep, by="Proteins", all.x = T)
    final_result <- almost[,c(1,dim(almost)[2],2:(dim(almost)[2]-1)) ]
    
    #OUTPUT: spitting results
    print("OUTPUT: list with three objects: protein list, peptide list, removed contaminants list")
    return(list(protBiorepIntensity, peptideBiorepIntensity, contaminants))
    # protBiorepIntensity <<- protBiorepIntensity
    # peptideBiorepIntensity <<- peptideBiorepIntensity
    # contaminants <<- contaminants
    # final_result <<- final_result
    # proteinGroups <<- proteinGroups
    # annotation <<- keys
    # evidence <<- data
}


create_saintq_files <- function(protBiorepIntensity, peptideBiorepIntensity, final_result, keys){
    #KEYS is from reading the annotation
      
    # SAINTQ HEADER (two extra rows on top)
    x <- t(keys[,c('SAINT','Condition','BioReplicaSaint')])
    extra <- cbind(c('','','Proteins'),c('','','Sequence'))
    header <- t(cbind(extra, x))
    theader <- t(header)
    checkthis <- data.frame(theader, row.names = NULL, stringsAsFactors = F)
    names(checkthis) = checkthis[3,]
    
    
    # PROTEINS: ADDING HEADER, merging based on row names
    proteinssaintqheader <- rbind(checkthis, final_result)
    
    # SEQUENCES: Adding HEADER:
    names(peptideBiorepIntensity)[2]<-"Sequence" #fixing the peptide feature name. to be changed earlier in the code maybe?
    sequencesaintqheader <- rbind(checkthis, peptideBiorepIntensity)
    
    #ADDED BY JUAN: to ensure replicates are in adjacent columns
    sequencesaintqheader<-cbind(sequencesaintqheader[, 1:2], 
                                sequencesaintqheader[, order(colnames(
                                    sequencesaintqheader[, 3:dim(sequencesaintqheader)[2]]))+2])
    
    proteinssaintqheader<-cbind(proteinssaintqheader[, 1:2], 
                                proteinssaintqheader[, order(colnames(
                                    proteinssaintqheader[, 3:dim(proteinssaintqheader)[2]]))+2])
    
    # Writing the files for PROTEIN
    output <- 'saintq/saintq_input_proteins.txt' # be careful changing this name because you also have to update the config file
    dir.create("saintq", showWarnings = F)
    write.table(proteinssaintqheader, output, sep='\t', row.names=F, col.names=F, quote=F)
    
    
    cat ("
         ### SAINTq parameter file
         ## use # to mark a line as comment
         
         ## normalize control intensities
         normalize_control=false
         
         ## name of file with intensities
         input_filename=saintq_input_proteins.txt
         
         ## valid: protein, peptide, fragment
         input_level=protein
         
         ## column names
         protein_colname=Proteins
         
         ## control bait selection rules
         compress_n_ctrl=100
         
         ## test bait replicate selection rules
         compress_n_rep=100
         ", file="saintq/saintq-config-proteins")
    
    
    # WRITING the files for SEQUENCES
    outsequences <- 'saintq/saintq_input_peptides.txt' # be careful changing this name because you also have to update the config file
    write.table(sequencesaintqheader, outsequences, sep = '\t', row.names = F, col.names = F, quote = F)
    
    cat ("
         ### SAINTq parameter file
         ## use # to mark a line as comment
         
         ## normalize control intensities
         normalize_control=false
         
         ## name of file with intensities
         input_filename=saintq_input_peptides.txt
         
         ## type of intensity
         ## valid: protein, peptide, fragment
         input_level=peptide
         
         ## column names
         protein_colname=Proteins
         pep_colname=Sequence
         
         
         ## control bait selection rules
         compress_n_ctrl=100
         
         ## test bait replicate selection rules
         compress_n_rep=100
         
         ## peptide selection rules
         min_n_pep=3
         best_prop_pep=0.5
         
         ", file="saintq/saintq-config-peptides")
    
    cat("\nDone! Check inside the new 'saintq' folder. \nYou should find 4 files:\n\n")
    cat("\t- saintq-config-peptides\n")
    cat("\t- saintq-config-proteins\n")
    cat("\t- saintq_input_peptides.txt\n")
    cat("\t- saintq_input_proteins.txt\n\n")
    cat("Now 'cd saintq' and run 'saintq saintq-config-peptides' and 'saintq saintq-config-proteins' and get the results right away\n\n")
    cat("You are welcome!\n\n")
}


# 
# ####Running some tests to compare scores####
# inputNewPeptide<-read.delim("saintq_new/saintq_input_peptides.txt")
# inputOldPeptide<-read.delim("saintq_old/saintq_input_peptides.txt")
# inputNewProtein<-read.delim("saintq_new/saintq_input_proteins.txt")
# inputOldProtein<-read.delim("saintq_old/saintq_input_proteins.txt")
# 
# colnames(inputNewPeptide)<-as.character(unlist(inputNewPeptide[2,])) #converting row names
# inputNewPeptide<-inputNewPeptide[-c(1,2),]
# 
# #a sample column comparison
# colNew<-inputNewProtein[3:nrow(inputNewProtein),c(1,3)]
# colOld<-inputOldProtein[3:nrow(inputOldProtein),c(1,3)]
# m<-merge(colNew, colOld, by="X", suffixes = c("_new", "_old"))
# m$T_new<-as.numeric(as.character(m$T_new))
# m$T_old<-as.numeric(as.character(m$T_old))
