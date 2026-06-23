#################### ISOFORM SWITCHING ANALYSIS ###########################

########### Performed on quant.sf files from Salmon ########################

######### Following IsoformSwitchAnalyzeR tutorial ##########################

######### Installing required packages #################

BiocManager::install("IsoformSwitchAnalyzeR")

########### Load the package #############

library(IsoformSwitchAnalyzeR)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)


packageVersion('IsoformSwitchAnalyzeR') ### Check the package version

getwd()

#_______________________________________________________________________________

##################### STEP-1: IMPORT THE DATA #############################

########## Import Salmon data ################

salmonQuant <- importIsoformExpression(
  parentDir = "D:/Projects/Isoform_Switching/salmon_output/quant_files",
  addIsofomIdAsColumn = TRUE)


#################### Make Design Matrix ###########################

############## Load the metadata info ############################

metadata <- read.csv("D:/Projects/Isoform_Switching/Updated_Metadata.csv")

View(metadata)

##################### Create Design Matrix ################################

myDesign <- data.frame(
  sampleID  = colnames(salmonQuant$abundance)[-1],
  condition = metadata$condition[match(colnames(salmonQuant$abundance)[-1], metadata$sample)],
  patient   = metadata$patient_id[match(colnames(salmonQuant$abundance)[-1], metadata$sample)] 
)

######### Use the match function to assign correct conditions #############

print(myDesign)


##################### Create SwitchAnalyzeR List ############################

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "D:/Projects/Isoform_Switching/salmon_quant/ref_files/gencode.v44.annotation.gtf.gz",
  isoformNtFasta       = "D:/Projects/Isoform_Switching/salmon_quant/ref_files/gencode.v44.transcripts.fa.gz",
  showProgress = FALSE)


head(aSwitchList$isoformFeatures,2)

head(aSwitchList$exons,2)

head(aSwitchList$ntSequence,2)


#_______________________________________________________________________________

######################### STEP-2: PREFILTERING ############################

######### Remove irrelevant genes/isoforms or non-expressed isoforms #######

?preFilter

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE,
  reduceFurtherToGenesWithConsequencePotential = TRUE)

#_______________________________________________________________________________

#################### PERFORM PCA FOR INITIAL QUALITY ANALYSIS ################

################### Extract Abundance (TPM) ##########################

pca_counts <- as.matrix(SwitchListFiltered$isoformRepExpression[, -1])

pca_counts <- apply(pca_counts, 2, as.numeric)

rownames(pca_counts) <- SwitchListFiltered$isoformRepExpression$isoform_id

################## Filter: top 25% most variable isoforms #########################

row_vars <- apply(pca_counts, 1, var)

pca_counts_filtered <- pca_counts[row_vars > quantile(row_vars, 0.75), ]

################## Log Transform #########################

log_data <- t(log2(pca_counts_filtered + 1))

######################### Run PCA ###########################

pca_analysis <- prcomp(log_data, scale. = FALSE, center = TRUE)

############## Prepare data for plot ###################

pca_df <- as.data.frame(pca_analysis$x)

pca_df$SampleID <- rownames(pca_df)

######### Matching sample condition & sampleID ###################

pca_df$condition <- SwitchListFiltered$designMatrix$condition[
  match(pca_df$SampleID, SwitchListFiltered$designMatrix$sampleID)]

pca_df$patient <- SwitchListFiltered$designMatrix$patient[
  match(pca_df$SampleID, SwitchListFiltered$designMatrix$sampleID)]


################ PCA Plot #####################

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  # Keep the points completely standard and clean
  geom_point(size = 4) +    
  
  # Crucial: Use the patient column (HC001, HC002) for the text labels instead of SampleID
  geom_text_repel(aes(label = patient), size = 2, fontface = "bold", color = "gray20") +  
  
  # Your original high-contrast colors
  scale_color_manual(values = c("normal" = "dodgerblue", "tumor" = "firebrick")) +
  
  labs(title = "PCA of Isoform Expression (TPM)",
       x = sprintf("PC1 (%1.1f%%)", summary(pca_analysis)$importance[2,1] * 100),
       y = sprintf("PC2 (%1.1f%%)", summary(pca_analysis)$importance[2,2] * 100),
       color = "Condition") + 
  theme_bw()


########### Scree plot — how much variance each PC captures ###################

scree_df <- data.frame(PC = 1:min(10, length(pca_analysis$sdev)),
  Variance = summary(pca_analysis)$importance[2, 1:min(10, length(pca_analysis$sdev))] * 100)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line() + geom_point() +
  labs(title = "Scree Plot", y = "% Variance Explained", x = "Principal Component") +
  theme_bw()

##################### PLOT TPM ##############################

############### Reshape to long format ####################

tpm_long <- pca_counts_filtered %>%
  as.data.frame() %>%
  tibble::rownames_to_column("isoform_id") %>%
  pivot_longer(cols = -isoform_id, names_to = "SampleID", values_to = "TPM")

tpm_long <- tpm_long %>%
  left_join(SwitchListFiltered$designMatrix[, c("sampleID", "condition")],
    by = c("SampleID" = "sampleID"))

ggplot(tpm_long, aes(x = reorder(SampleID, condition), y = log2(TPM + 1), fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("normal" = "dodgerblue", "tumor" = "firebrick")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Log2(TPM + 1)", 
       title = "TPM Distribution per Sample", fill = "Condition")


#_______________________________________________________________________________

######################## STEP-3: IDENTIFY ISOFORM SWITCHES ##################

?isoformSwitchTestDEXSeq

SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE)


extractSwitchSummary(SwitchListAnalyzed)

#### 420 genes, 676 transcripts/isoforms & 464 switches 

summary(SwitchListAnalyzed)

############### Extract the genes with the biggest switches #################

top_switches <- extractTopSwitches(SwitchListAnalyzed, n = 10)

View(top_switches)

################## Plot the switches for top genes #######################

switchPlot(SwitchListAnalyzed, 
           gene = top_switches$gene_id[1])

#_______________________________________________________________________________

##################### STEP-4: ANALYZING ORFs ###################################

################## Only Known Isoforms ##########################

"orfAnalysis" %in% names(SwitchListAnalyzed)

##### Returns TRUE, meaning CDS from GTF have been applied already #######

#_______________________________________________________________________________

############# STEP-5: EXTRACT NUCLEOTIDE & AMINO ACID SEQUENCES ################

?extractSequence

SwitchListAnalyzed <- extractSequence(SwitchListAnalyzed, 
                                      removeLongAAseq=TRUE,
                                      alsoSplitFastaFile=TRUE,
                                      pathToOutput = "switchanalyzer_output")

summary(SwitchListAnalyzed)

#_______________________________________________________________________________

################ STEP-6: RUN EXTERNAL SEQUENCE ANALYSIS TOOLS ##################

########### Using fasta files obtained in extractsequence step ##############

###### CPC2 (webserver) #######

###### SIGNALP (webserver) ##########

###### NETSURFP-3 (webserver) ########

##### Alternatively, you can use IUPred2A (webserver) ######

##### Pfam DONE (Done through galaxy.eu PfamScan tool) #######

##### DeepLoc2 DONE (Webserver) ##################

##### DeepTMHMM DONE (webserver) ##############


#_______________________________________________________________________________

################# STEP-7: IMPORTING EXTERNAL SEQUENCE ANALYSIS ###############

############### Add CPC2 Analysis ####################

SwitchListAnalyzed <- analyzeCPC2(
  switchAnalyzeRlist   = SwitchListAnalyzed,
  pathToCPC2resultFile = "switchanalyzer_output/result_cpc2.txt",
  removeNoncodinORFs   = TRUE)


summary(SwitchListAnalyzed)


################ Add Pfam Analysis #########################

SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = SwitchListAnalyzed,
  pathToPFAMresultFile = "switchanalyzer_output/pfam_result.txt",
  showProgress=TRUE)

head(SwitchListAnalyzed$isoformFeatures$isoform_id)



############### Add IUPred2A Analysis ########################

SwitchListAnalyzed <- analyzeIUPred2A(
  switchAnalyzeRlist        = SwitchListAnalyzed,
  pathToIUPred2AresultFile = "switchanalyzer_output/IUPred2A_result.txt",
  showProgress = TRUE)

SwitchListAnalyzed


############### Add SignalP Analysis ########################

SwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = SwitchListAnalyzed,
  pathToSignalPresultFile  = "switchanalyzer_output/SignalP_result.txt")


################# Add DeepLoc2 Analysis ############################

deeploc_files <- c(
  "switchanalyzer_output/DeepLoc_1.csv", 
  "switchanalyzer_output/DeepLoc_2.csv",
  "switchanalyzer_output/DeepLoc_3.csv", 
  "switchanalyzer_output/DeepLoc_4.csv")

SwitchListAnalyzed <- analyzeDeepLoc2(
  switchAnalyzeRlist = SwitchListAnalyzed,
  pathToDeepLoc2resultFile = deeploc_files,
  quiet = FALSE)

################### Add DeepTMHMM Analysis #########################

deeptmhmm_files <- c(
  "switchanalyzer_output/DeepTMHMM_1.gff3", 
  "switchanalyzer_output/DeepTMHMM_2.gff3",
  "switchanalyzer_output/DeepTMHMM_3.gff3",
  "switchanalyzer_output/DeepTMHMM_4.gff3")

SwitchListAnalyzed <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = SwitchListAnalyzed,
  pathToDeepTMHMMresultFile = deeptmhmm_files,
  showProgress=TRUE)

#_______________________________________________________________________________

################# STEP-8: PREDICT ALTERNATIVE SPLICING ###################

SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchListAnalyzed,
  quiet=FALSE)

table(SwitchListAnalyzed$AlternativeSplicingAnalysis$IR)

#_______________________________________________________________________________

####################### STEP-9: PREDICT SWITCH CONSEQUENCES ####################

?analyzeSwitchConsequences
SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze = "all", 
  dIFcutoff = 0.05,
  showProgress=TRUE)


extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = FALSE)

extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = TRUE)

extractConsequenceSummary(SwitchListAnalyzed)

summary(abs(SwitchListAnalyzed$isoformSwitchAnalysis$dIF))

SwitchListAnalyzed

table(SwitchListAnalyzed$isoformFeatures$switchConsequence)

#_______________________________________________________________________________

################### STEP-10: POST ANALYSIS ########################

######### Analyzing individual isoform switching ############

##### Genes/Isoforms with largest change in isoform usage #########

#### Either by smallest q-values (genes with most significant switches) #######

#### Or by largest dIF values (genes with largest effect size) ##########

######## Extract top switching genes by q-value #############

qvalue_cutoff <- 0.05

dif_cutoff <- 0.1

top_switch_conseq_q <- extractTopSwitches(
  SwitchListAnalyzed,
  filterForConsequences = TRUE,
  n = 10,
  sortByQvals = TRUE, 
  dIFcutoff = dif_cutoff,
  alpha = qvalue_cutoff) 

top_switch_conseq_q

write.csv(top_switch_conseq_q, "results/top_switch_conseq_q.csv")


########## Extract top switching genes by dif-value #############

top_switch_conseq_dif <- extractTopSwitches(
  SwitchListAnalyzed,
  filterForConsequences = TRUE,
  n = 10,
  sortByQvals = FALSE,
  dIFcutoff = dif_cutoff,
  alpha = qvalue_cutoff)

top_switch_conseq_dif

write.csv(top_switch_conseq_dif, "results/top_switch_conseq_dif.csv")


######## Extract data frame with all switching isoforms ###############

### Extract data.frame with all significant switches isoforms


allSwitches <- SwitchListAnalyzed$isoformFeatures

allSwitches

significantSwitches <- allSwitches[
  !is.na(allSwitches$isoform_switch_q_value) &
    allSwitches$isoform_switch_q_value < qvalue_cutoff &
    abs(allSwitches$dIF) > dif_cutoff, ]

nrow(significantSwitches)

length(unique(significantSwitches$gene_name))

write.csv(significantSwitches, "results/significant_switches_filtered.csv", row.names = FALSE)

View(significantSwitches)

#_______________________________________________________________________________

######################## STEP-11: VISUALIZATION ###############################

########## Volcano Like Plot ###############

ggplot(data = SwitchListAnalyzed$isoformFeatures,
       aes(x = gene_log2_fold_change, y = dIF)) +
  geom_point(
    aes(color = abs(dIF) > dif_cutoff & isoform_switch_q_value < qvalue_cutoff),
    size = 1) +
  facet_wrap(~ condition_2) +
  geom_hline(yintercept =  dif_cutoff, linetype = 'dashed', color = 'gray50') +
  geom_hline(yintercept = -dif_cutoff, linetype = 'dashed', color = 'gray50') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual('Significant\nIsoform Switch', values = c('black', 'red')) +
  labs(x = 'Gene log2 fold change', y = 'dIF') +
  theme_bw()



############ Switch Plot for the most significant gene ################

switchPlot(SwitchListAnalyzed, gene = 'PRX')

switchPlot(SwitchListAnalyzed, gene = "LAMA2")

switchPlot(SwitchListAnalyzed, gene = "MAD2L2")

switchPlot(SwitchListAnalyzed, gene = "FBLN5")

switchPlot(SwitchListAnalyzed, gene = "FBLN2")


########## Plot summarizing Alternative Splicing events ##################

extractSplicingSummary(SwitchListAnalyzed,
                       asFractionTotal = FALSE,
                       plotGenes=FALSE)

############ Splicing Enrichment graph ###################

splicingEnrichment <- extractSplicingEnrichment(SwitchListAnalyzed,
                                                splicingToAnalyze='all',
                                                returnResult=TRUE,
                                                returnSummary=TRUE)


extractSplicingGenomeWide(
  SwitchListAnalyzed,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = c('A3','MES','ATSS'), 
  plot=TRUE,
  returnResult=FALSE  
)



################ Analyze Biological Mechanisms behind Isoform Siwtching #############

##### Analyze the biological mechanisms ###########

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  SwitchListAnalyzed, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here


####### Subset to those with differences #########

bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]


### Extract the consequences of interest already stored in the switchAnalyzeRlist

myConsequences <- SwitchListAnalyzed$switchConsequence

myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]

myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated) # id for specific iso comparison


### Obtain the mechanisms of the isoform switches with consequences

bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)

bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  # id for specific iso comparison



### Create list with the isoPair ids for each consequence

AS   <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_structure')]
aTSS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tss'             )]
aTTS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tts'             )]

mechList <- list(
  AS=AS,
  aTSS=aTSS,
  aTTS=aTTS
)

### Create Venn diagram
library(VennDiagram)

myVenn <- venn.diagram(
  x = mechList,
  col='transparent',
  alpha=0.4,
  fill=RColorBrewer::brewer.pal(n=3,name='Dark2'),
  filename=NULL
)

### Plot the venn diagram
grid.newpage() ; grid.draw(myVenn)





