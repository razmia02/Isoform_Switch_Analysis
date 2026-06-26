################## INITIAL QC AFTER SALMON QUANTIFICATION ##################


################# Calling Required packages and libraries ################

library(dplyr)
library(tidyr)
library(ggplot2)

################## Read in the data ######################

mapping_df <- read.csv("D:/Projects/Isoform_Switching/Mapping_Rate.csv")

####### Check what's in the column ##########

str(mapping_df)
head(mapping_df$mapping_rate)

########## Convert to numeric ###############

mapping_df$mapping_rate <- as.numeric(gsub("%", "", mapping_df$mapping_rate))

##################### Plot Mapping rate ######################

ggplot(mapping_df, aes(x = reorder(sample, mapping_rate), 
                       y = mapping_rate, 
                       fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 75, linetype = "dashed", 
             color = "red", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", mapping_rate)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  coord_flip() +
  labs(
    title   = "Salmon Mapping Rate per Sample",
    x       = NULL,
    y       = "Mapping Rate (%)",
    fill    = "Condition",
    caption = "Red dashed line = 75% minimum threshold"
  ) +
  theme_bw() +
  theme(plot.title  = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))


################ Read in the quant files data ######################

salmon_dir <- "D:/Projects/Isoform_Switching/salmon_output/quant_files"

########## Find all quant.sf files ###########

quant_files <- list.files(salmon_dir, 
                          pattern = "quant.sf", 
                          recursive = TRUE, 
                          full.names = TRUE)

########## Extract total counts per sample #################

count_summary <- lapply(quant_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t")
  sample_name <- basename(dirname(f))
  data.frame(
    sample        = sample_name,
    total_counts  = sum(df$NumReads),
    num_transcripts_detected = sum(df$NumReads > 0)
  )
}) %>% bind_rows()

print(count_summary)



########### Add condition from metadata ################

metadata <- read.csv("D:/Projects/Isoform_Switching/Paper/Supplementary_Table_1.csv")

print(metadata)


count_summary$condition <- metadata$condition[
  match(count_summary$sample, metadata$sample)]

######## Plot Total Mapped Counts per Sample #################

ggplot(count_summary, aes(x = reorder(sample, total_counts), 
                          y = total_counts / 1e6,  # convert to millions
                          fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1fM", total_counts / 1e6)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  scale_y_continuous(limits = c(0, max(count_summary$total_counts/1e6) * 1.15),
                     expand = c(0, 0)) +
  coord_flip() +
  labs(
    title   = "Total Mapped Read Counts per Sample",
    x       = NULL,
    y       = "Total Mapped Reads (Millions)",
    fill    = "Condition"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))


########## Plot Number of Transcripts Detected per Sample ##############

ggplot(count_summary, aes(x = reorder(sample, num_transcripts_detected), 
                          y = num_transcripts_detected / 1000,  # convert to thousands
                          fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1fK", num_transcripts_detected / 1000)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  scale_y_continuous(limits = c(0, max(count_summary$num_transcripts_detected/1000) * 1.15),
                     expand = c(0, 0)) +
  coord_flip() +
  labs(
    title   = "Number of Transcripts Detected per Sample",
    x       = NULL,
    y       = "Transcripts Detected (Thousands)",
    fill    = "Condition"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

ggsave("results/transcripts_detected.png", width = 8, height = 6, dpi = 300)

###### Plot Both metrics combined in one figure #############


count_long <- count_summary %>%
  mutate(
    `Mapped Reads (M)`        = total_counts / 1e6,
    `Transcripts Detected (K)` = num_transcripts_detected / 1000
  ) %>%
  select(sample, condition, `Mapped Reads (M)`, `Transcripts Detected (K)`) %>%
  pivot_longer(cols = c(`Mapped Reads (M)`, `Transcripts Detected (K)`),
               names_to = "metric", values_to = "value")

ggplot(count_long, aes(x = reorder(sample, value), 
                       y = value, 
                       fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  facet_wrap(~ metric, scales = "free_x") +  # free scales since units differ
  coord_flip() +
  labs(
    title = "Salmon Quantification QC Metrics",
    x     = NULL,
    y     = "Value",
    fill  = "Condition"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )



########### Merge mapping rates with count_summary ##################
combined_qc <- count_summary %>%
  left_join(mapping_df, by = c("sample", "condition"))



combined_long <- combined_qc %>%
  mutate(
    `Mapping Rate (%)` = mapping_rate,
    `Mapped Reads (M)` = total_counts / 1e6
  ) %>%
  select(sample, condition, `Mapping Rate (%)`, `Mapped Reads (M)`) %>%
  pivot_longer(cols = c(`Mapping Rate (%)`, `Mapped Reads (M)`),
               names_to = "metric", values_to = "value")

ggplot(combined_long, aes(x = reorder(sample, value), 
                          y = value, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  
  # Red threshold line only on mapping rate panel
  geom_hline(data = subset(combined_long, metric == "Mapping Rate (%)"),
             aes(yintercept = 75), 
             linetype = "dashed", color = "red", linewidth = 0.7) +
  
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  facet_wrap(~ metric, scales = "free_x") +
  coord_flip() +
  labs(
    title   = "Salmon Quantification QC",
    x       = NULL,
    y       = NULL,
    fill    = "Condition",
    caption = "Dashed line = 75% mapping rate threshold"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )


############ Merge both dataframes #############3

combined_qc <- count_summary %>%
  left_join(mapping_df, by = c("sample", "condition")) %>%
  mutate(
    `Mapping Rate (%)` = mapping_rate,
    `Mapped Reads (M)` = total_counts / 1e6,
    `Transcripts Detected (K)` = num_transcripts_detected / 1000
  ) %>%
  select(sample, condition, 
         `Mapping Rate (%)`, 
         `Mapped Reads (M)`, 
         `Transcripts Detected (K)`)

######### Pivot to long format ############

combined_long <- combined_qc %>%
  pivot_longer(
    cols      = c(`Mapping Rate (%)`, `Mapped Reads (M)`, `Transcripts Detected (K)`),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  # Control panel order
  mutate(metric = factor(metric, levels = c("Mapping Rate (%)", 
                                            "Mapped Reads (M)", 
                                            "Transcripts Detected (K)")))

write.csv(combined_long, "results/combined_qc_metrics.csv")

############# Plot all three metrics combined ################

ggplot(combined_long, aes(x = reorder(sample, value), 
                          y = value, 
                          fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  
  # Threshold line ONLY on mapping rate panel
  geom_hline(
    data      = subset(combined_long, metric == "Mapping Rate (%)"),
    aes(yintercept = 75),
    linetype  = "dashed", 
    color     = "red", 
    linewidth = 0.7
  ) +
  
  # Value labels at end of each bar
  geom_text(aes(label = ifelse(metric == "Mapping Rate (%)",
                               sprintf("%.1f%%", value),
                               sprintf("%.1f",   value))),
            hjust  = -0.1, 
            size   = 2.5) +
  
  scale_fill_manual(values = c("normal" = "dodgerblue", 
                               "tumor"  = "firebrick")) +
  
  # Free scales — each panel has different units
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  
  coord_flip() +
  
  labs(
    title   = "Salmon Quantification QC Metrics",
    x       = NULL,
    y       = NULL,
    fill    = "Condition",
    caption = "Dashed line = 75% mapping rate threshold"
  ) +
  
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", size = 13),
    strip.text   = element_text(face = "bold", size = 10),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    legend.position = "bottom",
    panel.spacing   = unit(1, "lines"),  # space between panels
    plot.caption    = element_text(color = "gray50")
  )


