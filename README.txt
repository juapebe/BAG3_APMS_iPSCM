###NOTES###

##Experimental design
Four replicates. 14 different conditions. Each one comes from a different cell line.
The value we are about is "Intensity". Number of peptide counts can also be of help.
For cardiomyocytes, every cell line was differentiated multiple times, and pooled together to give enough for an experiment. Cells were frozen after differentiation.
The actual affinity purification step (most variability introduced here) was done in parallel for each one of the replicates (14 samples per day).
The mass spectrometry was done for all at the same time.
Quality control: 
+BAG3 should not happen in the controls, but it does show in at least one of them. Intensity of BAG3 should be a proxy for how well the pulldown worked (i.e. all the other interactions should be dependant on BAG3 intensity levels). 
+ 11 peptides were spiked in. These should be a control for how well the mass spec performed. All labeled/aggregated into fictional protein "Biognosys|iRT-Kit_WR_fusion".
Most common charge for these is 2, so this should probably the charge used for normalization (sequences ending in "__2")

##Aims.
Multiple questions can be asked
Start with *'compare interactions of BAG3 between 'wild-type' (healthy) and disease allele (K22),in cardiomyocytes*. (that is, A4 vs K22)
Also, what are legitimate interactors of BAG3 (that is, WTc11 vs A4) 
Not definite analysis - just come up with a list of interactors that then can be finely quantified using targeted proteomics.

##Files in folder
experimentDescription: a small diagram on how everything was done.
finalAPMS_dataAnalysis: general handling of data
auxiliaryFunctions.R: side functions to perform small tasks.
process4MSStats_3.R: for processing files for MSStats
process4saintq_2.R: for processing files for SAINTQ (using parse_maxquant to load files)
SetA+B+C+D_Fusion: data.
saintq: saintq files
MSStats: MSStats files

##Conditions
As in "annotation.txt" file:
+WTc11 - control cell line. Cardiomyocytes
+A4 - wild-type BAG3. Cardiomyocytes
+C10 - 'protective' BAG3 allele. Cardiomyocytes
+P9 - 'really bad' BAG3 allele, expected gain-of-function. Cardiomyocytes
+K22 - 'disease' BAG3 allele. Cardiomyocytes
+homo+DMSO - wild-type BAG3 (both alleles) treated with control drug. Cardiomyocytes
+homo+Bort - wild-type BAG3 (both alleles) treated with test drug. Cardiomyocytes
+HSPB7 - HSPB7 pull down (not BAG3!)
+TetOn-BAG3 - BAG3 massively overexpressed. Cardiomyocytes
+iPS_WTc11 - control cell line. iPS cells
+iPS_A4 - wild-type BAG3. iPS cells
+iPS_C10 - 'protective' BAG3 allele. iPS cells
+iPS_P9 - 'really bad' BAG3 allele, expected gain-of-function. iPS cells
+iPS_K22 - 'disease' BAG3 allele. iPS cells
+RBM20 (only one replicate. test. ignore)