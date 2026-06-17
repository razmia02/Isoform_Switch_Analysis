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
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = metadata$condition[match(colnames(salmonQuant$abundance)[-1]
                                       , metadata$sample)])

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
  removeSingleIsoformGenes = TRUE, 
  reduceFurtherToGenesWithConsequencePotential = TRUE)

#_______________________________________________________________________________

#################### PERFORM PCA FOR INITIAL QUALITY ANALYSIS ################

################### Extract Abundance (TPM) ##########################

pca_counts <- SwitchListFiltered$isoformRepExpression[, -1]

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

pca_df$Condition <- SwitchListFiltered$designMatrix$condition[
  match(pca_df$SampleID, SwitchListFiltered$designMatrix$sampleID)]

################ PCA Plot #####################

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = SampleID)) +
  geom_point(size = 4) +    
  geom_text_repel(size = 3) +  
  scale_color_manual(values = c("dodgerblue", "firebrick")) +
  labs(title = "PCA of Isoform Expression (TPM)",
    x = sprintf("PC1 (%1.1f%%)", summary(pca_analysis)$importance[2,1] * 100),
    y = sprintf("PC2 (%1.1f%%)", summary(pca_analysis)$importance[2,2] * 100),
    color = "Condition") +
  theme_bw() +
  theme(legend.position = "right")



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

tpm_long <- pca_counts %>%
  tibble::rownames_to_column("isoform_id") %>%
  pivot_longer(cols = -isoform_id, names_to = "SampleID", values_to = "TPM")

tpm_long <- tpm_long %>%
  left_join(SwitchListFiltered$designMatrix[, c("sampleID", "condition")],
    by = c("SampleID" = "sampleID"))

ggplot(tpm_long, aes(x = SampleID, y = log2(TPM + 1), fill = condition)) +
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

#### 475 genes, 773 transcripts/isoforms & 542 switches 

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

"orfAnalysis" %in% names(aSwitchList)

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
  pathToPFAMresultFile = "switchanalyzer_output/pfam_results_galaxy.txt",
  showProgress=TRUE)

head(SwitchListAnalyzed$isoformFeatures$isoform_id)



############### Add IUPred2A Analysis ########################

SwitchListAnalyzed <- analyzeIUPred2A(
  switchAnalyzeRlist        = SwitchListAnalyzed,
  pathToIUPred2AresultFile = "switchanalyzer_output/results_IUPred2A.txt",
  showProgress = TRUE)

SwitchListAnalyzed


############### Add SignalP Analysis ########################

SwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = SwitchListAnalyzed,
  pathToSignalPresultFile  = "switchanalyzer_output/signalp_results.txt")

#_______________________________________________________________________________

################# STEP-8: PREDICT ALTERNATIVE SPLICING ###################

SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchListAnalyzed,
  quiet=FALSE)

table(SwitchListAnalyzed$AlternativeSplicingAnalysis$IR)

#_______________________________________________________________________________

####################### STEP-9: PREDICT SWITCH CONSEQUENCES ####################

consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.1,
  showProgress=TRUE)


extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = FALSE)

extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = TRUE)

extractConsequenceSummary(SwitchListAnalyzed)

#_______________________________________________________________________________

################### STEP-10: POST ANALYSIS ########################

######### Analyzing individual isoform switching ############

##### Genes/Isoforms with largest change in isoform usage #########

#### Either by smallest q-values (genes with most significant switches) #######

#### Or by largest dIF values (genes with largest effect size) ##########

######## Extract top switching genes by q-value #############

top_switch_conseq_q <- extractTopSwitches(SwitchListAnalyzed, 
                   filterForConsequences = TRUE, n = 10, 
                   sortByQvals = TRUE)

top_switch_conseq_q

write.csv(top_switch_conseq_q, "results/top_switch_conseq_q.csv")


########## Extract top switching genes by q-value #############

top_switch_conseq_dif <- extractTopSwitches(SwitchListAnalyzed, 
                                            filterForConsequences = TRUE, 
                                            n = 10, 
                                            sortByQvals = FALSE)

write.csv(top_switch_conseq_dif, "results/top_switch_conseq_dif.csv")


######## Extract data frame with all switching isoforms ###############

### Extract data.frame with all switching isoforms
switchingIso <- extractTopSwitches(SwitchListAnalyzed, 
                                   filterForConsequences = TRUE, 
                                   n = NA,
                                   extractGenes = TRUE,  
                                   sortByQvals = TRUE)

nrow(switchingIso) ### 437 genes have functional consequences switching

write.csv(switchingIso, "results/All_Switches.csv")

subset(switchingIso, gene_name == 'CRTC1')


#_______________________________________________________________________________

######################## STEP-11: VISUALIZATION ###############################

########## Volcano Like Plot ###############

ggplot(data=SwitchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()


############ Switch Plot for DIXDC1 gene ######################

switchPlot(SwitchListAnalyzed, gene = 'CRTC1')

switchPlotIsoUsage(SwitchListAnalyzed, gene = 'CRTC1')

########## Plot summarizing Alternative Splicing events ##################

extractSplicingSummary(SwitchListAnalyzed,
                       asFractionTotal = FALSE,
                       plotGenes=FALSE)

############ Splicing Enrichment graph ###################

splicingEnrichment <- extractSplicingEnrichment(SwitchListAnalyzed,
                                                splicingToAnalyze='all',
                                                returnResult=TRUE,
                                                returnSummary=TRUE)










