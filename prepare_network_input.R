###################################################
# title: Prepare network input                    #
# author: Marie-Madlen Pust                       #
# date: 15 12 2020                                #
##################################################


########################################################
#clear environment
rm(list=ls())

#load packages
library('readr')
library('dplyr')
library('plyr')
library('tidyr')
library('Hmisc')
library('reshape2')
library('stringr')
library('matrixStats')

# most abundant species
abund_spec = 99.9
weight_val = 0.60
weight_val2 = -1.10
sig_level = 0.01

########################################################
# cough swabs (healthy) ####
# import data files
df_throat_complete <- read_delim(
  "input_files/illumina_healthy_controls.csv", 
  ";", escape_double = FALSE, trim_ws = TRUE)
df_throat_complete <- data.frame(df_throat_complete)
df_throat_complete$Genus <- NULL

df_throat_complete <- df_throat_complete[-1,]
df_throat_complete <- df_throat_complete[
  complete.cases(df_throat_complete), ]
df_throat_complete$Species <- str_replace(
  df_throat_complete$Species, " ", "_")
df_throat_complete <- ddply(
  df_throat_complete,"Species",numcolwise(sum))
df_throat_complete <- df_throat_complete[
  - grep("phage", df_throat_complete$Species),]
df_throat_complete <- df_throat_complete[
  - grep("virus", df_throat_complete$Species),]
df_throat_complete <- df_throat_complete[
  - grep("albicans", df_throat_complete$Species),]
df_throat_complete <- df_throat_complete[
  - grep("fumigatus", df_throat_complete$Species),]
df_throat_complete <- df_throat_complete[
  - grep("dubliniensis", df_throat_complete$Species),]
rownames(df_throat_complete) <- df_throat_complete$Species
df_throat_complete$Species <- NULL
df_throat_complete <- df_throat_complete[-1,]
percent_throat_complete <- round(
  (ncol(df_throat_complete) * 33.3) / 100, 0)
df_throat_complete <- df_throat_complete[
  rowSums(df_throat_complete == 0) <= percent_throat_complete, ]

water_blank <- c("MBCF_83", "MBCF_84", "MBCF_85", 
                 "MBCF_86", "MBCF_87", "MBCF_88")
swabs_blank <- c("MBCF_80", "MBCF_81", "MBCF_82", 
                 "MBCF_103")
df_controls <- c(water_blank, swabs_blank)

df_water_blank <- select(
  df_throat_complete, 
  water_blank)
df_swabs_blank <- select(
  df_throat_complete, 
  swabs_blank)  
df_healthy_counts <- select(
  df_throat_complete, 
  -df_controls)

# extract most abundant species
# healthy
sum_healthy = colSums(
  df_healthy_counts)
sum_all_healthy = sum(
  sum_healthy)
df_healthy_counts$abundance <- (rowSums(
  df_healthy_counts[,1:ncol(df_healthy_counts)]) / sum_all_healthy) * 100
# sort abundance decreasing
df_healthy_counts <- df_healthy_counts[
  with(df_healthy_counts, order(-abundance)), ]

df_healthy_counts$cumsum <- cumsum(
  df_healthy_counts$abundance)
# table with abundant species
df_healthy_counts <- subset(
  df_healthy_counts, 
  cumsum <= abund_spec)
df_healthy_counts$cumsum <- NULL
df_healthy_counts$abundance <- NULL

# extract absolute abundance values
df_healthy_abundance_value <- rowMedians(
  as.matrix(df_healthy_counts)) 
df_healthy_abundance_species <- rownames(df_healthy_counts)
df_healthy_node_information <- data.frame(
  cbind(df_healthy_abundance_value,
        df_healthy_abundance_species))
df_healthy_node_information <- separate(
  df_healthy_node_information,
  col = 'df_healthy_abundance_species',
  into=c('Genus', 'N'), sep='_', remove=TRUE)
df_healthy_node_information$N <- NULL

# convert columns from factor to numeric
df_healthy_node_information$df_healthy_abundance_value <-
  as.numeric(as.character(
    df_healthy_node_information$df_healthy_abundance_value))
df_healthy_node_information <- ddply(
  df_healthy_node_information,"Genus",numcolwise(sum))
df_healthy_node_information$df_healthy_abundance_value <- round(
  df_healthy_node_information$df_healthy_abundance_value,2)
colnames(df_healthy_node_information) <- c("Label", "Value")
df_healthy_counts <- data.frame(
  t(df_healthy_counts))

# generate correlation matrix
# healthy
df_healthy_counts_cor <- rcorr(
  as.matrix(df_healthy_counts), 
  type = 'spearman')

# create node and edge lists
# healthy
df_healthy_counts_cor_Pval <- df_healthy_counts_cor$P
df_healthy_counts_cor_COR <- df_healthy_counts_cor$r

# generate edgelist
df_healthy_counts_cor_edges <- melt(
  df_healthy_counts_cor_Pval)
df_healthy_counts_cor_coredges <- melt(
  df_healthy_counts_cor_COR)
df_healthy_counts_cor_edges$COR <- df_healthy_counts_cor_coredges$value
df_healthy_counts_cor_edges$value <- round(
  df_healthy_counts_cor_edges$value, 5)
df_healthy_counts_cor_edges$Label <- row.names(
  df_healthy_counts_cor_coredges)

colnames(df_healthy_counts_cor_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_healthy_counts_cor_edges) <- NULL
df_healthy_counts_cor_edges <- df_healthy_counts_cor_edges[
  complete.cases(df_healthy_counts_cor_edges), ]
df_healthy_counts_cor_edges_short <- subset(
  df_healthy_counts_cor_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_healthy_counts_cor_edges_short$cortype <- ifelse(
  df_healthy_counts_cor_edges_short$Weight > 0, "pos", "neg")

df_healthy_counts_cor_edges_short$Genus <- 
  df_healthy_counts_cor_edges_short$Source
df_healthy_counts_cor_edges_short <- separate(
  df_healthy_counts_cor_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_healthy_counts_cor_edges_short$N <- NULL
rownames(df_healthy_counts_cor_edges_short) <- NULL
df_healthy_counts_cor_edges_short$Weight <- round(
  df_healthy_counts_cor_edges_short$Weight,2)
df_healthy_counts_cor_edges_short$Type <- "Undirected"

df_healthy_counts_cor_nodes <- select(
  df_healthy_counts_cor_edges_short,
  c("Source", "Genus"))
df_healthy_counts_cor_nodes = df_healthy_counts_cor_nodes[
  !duplicated(df_healthy_counts_cor_nodes$Source),]
df_healthy_counts_cor_nodes <- data.frame(
  df_healthy_counts_cor_nodes)
rownames(df_healthy_counts_cor_nodes) <- NULL
colnames(df_healthy_counts_cor_nodes) <- c("Id", "Label")
df_healthy_counts_cor_nodes1 <- merge(
  df_healthy_counts_cor_nodes, 
  df_healthy_node_information[,c("Label", "Value")], 
  by="Label")

# extract most abundant species
# water controls
# keep species that are non-zero in half of the controls
percent_water_blank <- round(
  (ncol(df_water_blank) * 33.3) / 100, 0)
df_water_blank <- df_water_blank[
  rowSums(df_water_blank == 0) <= percent_water_blank, ]

sum_water = colSums(
  df_water_blank)
sum_all_water = sum(
  sum_water)
df_water_blank$abundance <- (rowSums(
  df_water_blank[,1:ncol(df_water_blank)]) / sum_all_water) * 100
# sort abundance decreasing
df_water_blank <- df_water_blank[
  with(df_water_blank, order(-abundance)), ]

df_water_blank$cumsum <- cumsum(
  df_water_blank$abundance)
# table with abundant species
df_water_blank <- subset(
  df_water_blank, 
  cumsum <= abund_spec)
df_water_blank$cumsum <- NULL
df_water_blank$abundance <- NULL
# extract absolute abundance values
df_water_abundance_value <- rowMedians(
  as.matrix(df_water_blank)) 
df_water_abundance_species <- rownames(
  df_water_blank)
df_water_node_information <- data.frame(
  cbind(df_water_abundance_value, 
        df_water_abundance_species))
df_water_node_information <- separate(
  df_water_node_information,
  col = 'df_water_abundance_species',
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_water_node_information$N <- NULL
# convert columns from factor to numeric
df_water_node_information$df_water_abundance_value <- as.numeric(
  as.character(df_water_node_information$df_water_abundance_value))
df_water_node_information <- ddply(
  df_water_node_information,
  "Genus",
  numcolwise(sum))
df_water_node_information$df_water_abundance_value <- round(
  df_water_node_information$df_water_abundance_value,2)
colnames(df_water_node_information) <- c("Label", "Value")
df_water_blank <- data.frame(
  t(df_water_blank))

# generate correlation matrix
# water
df_water_blank_cor <- rcorr(
  as.matrix(df_water_blank), type = 'spearman')

# Water
df_water_blank_cor_Pval <- df_water_blank_cor$P
df_water_blank_cor_COR <- df_water_blank_cor$r

# generate edgelist
df_water_blank_cor_edges <- melt(
  df_water_blank_cor_Pval)
df_water_blank_cor_edges_coredges <- melt(
  df_water_blank_cor_COR)
df_water_blank_cor_edges$COR <- df_water_blank_cor_edges_coredges$value
df_water_blank_cor_edges$value <- round(
  df_water_blank_cor_edges$value, 5)
df_water_blank_cor_edges$Label <- row.names(
  df_water_blank_cor_edges_coredges)

colnames(df_water_blank_cor_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_water_blank_cor_edges) <- NULL
df_water_blank_cor_edges <- df_water_blank_cor_edges[
  complete.cases(df_water_blank_cor_edges), ]
df_water_blank_cor_edges_short <- subset(
  df_water_blank_cor_edges, 
  pValue < sig_level & Weight >= weight_val)

df_water_blank_cor_edges_short$Genus <- df_water_blank_cor_edges_short$Source
df_water_blank_cor_edges_short <- separate(
  df_water_blank_cor_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_water_blank_cor_edges_short$N <- NULL
rownames(df_water_blank_cor_edges_short) <- NULL
df_water_blank_cor_edges_short$Weight <- round(
  df_water_blank_cor_edges_short$Weight,2)
df_water_blank_cor_edges_short$Type <- "Undirected"

df_water_blank_cor_nodes <- select(
  df_water_blank_cor_edges_short, 
  c("Source", "Genus"))
df_water_blank_cor_nodes = df_water_blank_cor_nodes[
  !duplicated(df_water_blank_cor_nodes$Source),]
df_water_blank_cor_nodes <- data.frame(
  df_water_blank_cor_nodes)
rownames(df_water_blank_cor_nodes) <- NULL
colnames(df_water_blank_cor_nodes) <- c("Id", "Label")
df_water_blank_cor_nodes1 <- merge(
  df_water_blank_cor_nodes, 
  df_water_node_information[, c("Label", "Value")], by="Label")

# extract most abundant species
# blank swabs
# keep species that are non-zero in half of the controls
percent_swab_blank <- round(
  (ncol(df_swabs_blank) * 33.3) / 100, 0)
df_swabs_blank <- df_swabs_blank[
  rowSums(df_swabs_blank == 0) <= percent_swab_blank, ]

sum_swab = colSums(
  df_swabs_blank)
sum_all_swabs = sum(
  sum_swab)

df_swabs_blank$abundance <- (
  rowSums(df_swabs_blank[,1:ncol(df_swabs_blank)]) / sum_all_swabs) * 100
# sort abundance decreasing
df_swabs_blank <- df_swabs_blank[
  with(df_swabs_blank, order(-abundance)), ]

df_swabs_blank$cumsum <- cumsum(
  df_swabs_blank$abundance)
# table with abundant species
df_swabs_blank <- subset(
  df_swabs_blank, cumsum <= abund_spec)
df_swabs_blank$cumsum <- NULL
df_swabs_blank$abundance <- NULL
# extract absolute abundance values
df_swabs_abundance_value <- rowMedians(
  as.matrix(df_swabs_blank)) 
df_swabs_abundance_species <- rownames(
  df_swabs_blank)
df_swabs_node_information <- data.frame(
  cbind(df_swabs_abundance_value,
        df_swabs_abundance_species))
df_swabs_node_information <- separate(
  df_swabs_node_information,
  col = 'df_swabs_abundance_species',
  into=c('Genus', 'N'),
  sep='_', 
  remove=TRUE)
df_swabs_node_information$N <- NULL
# convert columns from factor to numeric
df_swabs_node_information$df_swabs_abundance_value <- as.numeric(
  as.character(df_swabs_node_information$df_swabs_abundance_value))
df_swabs_node_information <- ddply(
  df_swabs_node_information,
  "Genus",
  numcolwise(sum))
df_swabs_node_information$df_swabs_abundance_value <- round(
  df_swabs_node_information$df_swabs_abundance_value,2)
colnames(df_swabs_node_information) <- c("Label", "Value")
df_swabs_blank1 <- data.frame(
  t(df_swabs_blank))
df_swabs_blank2 <- colMeans(df_swabs_blank1)
df_swabs_blank <- data.frame(
  rbind(df_swabs_blank1, df_swabs_blank2))

# Generate correlation matrix
# Blank
df_swabs_blank_cor <- rcorr(
  as.matrix(df_swabs_blank), 
  type = 'spearman')
df_swab_blank_cor_Pval <- df_swabs_blank_cor$P
df_swab_blank_cor_COR <- df_swabs_blank_cor$r

# generate edgelist
df_swab_blank_cor_edges <- melt(
  df_swab_blank_cor_Pval)
df_swab_blank_cor_edges_coredges <- melt(
  df_swab_blank_cor_COR)
df_swab_blank_cor_edges$COR <- df_swab_blank_cor_edges_coredges$value
df_swab_blank_cor_edges$value <- round(
  df_swab_blank_cor_edges$value, 5)
df_swab_blank_cor_edges$Label <- row.names(
  df_swab_blank_cor_edges_coredges)

colnames(df_swab_blank_cor_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_swab_blank_cor_edges) <- NULL
df_swab_blank_cor_edges <- df_swab_blank_cor_edges[
  complete.cases(df_swab_blank_cor_edges), ]
df_swab_blank_cor_edges_short <- subset(
  df_swab_blank_cor_edges, pValue < sig_level & Weight >= weight_val)

df_swab_blank_cor_edges_short$Genus <- df_swab_blank_cor_edges_short$Source
df_swab_blank_cor_edges_short <- separate(
  df_swab_blank_cor_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_swab_blank_cor_edges_short$N <- NULL
rownames(df_swab_blank_cor_edges_short) <- NULL
df_swab_blank_cor_edges_short$Weight <- round(
  df_swab_blank_cor_edges_short$Weight,2)
df_swab_blank_cor_edges_short$Type <- "Undirected"

df_swab_blank_cor_nodes <- select(
  df_swab_blank_cor_edges_short, c("Source", "Genus"))
df_swab_blank_cor_nodes = df_swab_blank_cor_nodes[
  !duplicated(df_swab_blank_cor_nodes$Source),]
df_swab_blank_cor_nodes <- data.frame(
  df_swab_blank_cor_nodes)
rownames(df_swab_blank_cor_nodes) <- NULL
colnames(df_swab_blank_cor_nodes) <- c("Id", "Label")
df_swab_blank_cor_nodes1 <- merge(
  df_swab_blank_cor_nodes, 
  df_water_node_information[,c("Label", "Value")], by="Label")

#################################################################
# SRR data ####
# import data files
df_SRR_complete <- read_delim(
  "input_files/illumina_SRR_sputum.csv", ";", 
  escape_double = FALSE, trim_ws = TRUE)
df_SRR_complete <- data.frame(
  df_SRR_complete)

df_SRR_complete <- ddply(
  df_SRR_complete,"Species", numcolwise(sum))
df_SRR_complete <- df_SRR_complete[-c(1:26),]
df_SRR_complete <- df_SRR_complete[
  - grep("phage", df_SRR_complete$Species),]
df_SRR_complete <- df_SRR_complete[
  - grep("virus", df_SRR_complete$Species),]
df_SRR_complete <- df_SRR_complete[
  - grep("Pichia", df_SRR_complete$Species),]
df_SRR_complete <- df_SRR_complete[
  - grep("albicans", df_SRR_complete$Species),]

df_SRR_complete <- df_SRR_complete[
  complete.cases(df_SRR_complete), ]
df_SRR_complete$Species <- str_replace(
  df_SRR_complete$Species, " ", "_")
rownames(df_SRR_complete) <- df_SRR_complete$Species
df_SRR_complete$Species <- NULL

SRR_complete_severe <- c(
  "SRR3284698_CF", "SRR3284703_CF", "SRR3286492_CF", "SRR3286493_CF", "SRR3286495_CF")
SRR_complete_mild <- c(
  "SRR3284701_CF", "SRR3284702_CF", "SRR3284706_CF", "SRR3286490_CF", "SRR3286491_CF", 
  "SRR3286494_CF")

df_SRR_complete_severe <- select(
  df_SRR_complete, SRR_complete_severe)

df_SRR_complete_mild <- select(
  df_SRR_complete, SRR_complete_mild)

# severe cases
sum_SRR_severe = colSums(
  df_SRR_complete_severe)
sum_all_SRR_severe = sum(
  sum_SRR_severe)
df_SRR_complete_severe$abundance <- (rowSums(
  df_SRR_complete_severe[,1:ncol(df_SRR_complete_severe)]) / sum_all_SRR_severe) * 100
# sort abundance decreasing
df_SRR_complete_severe <- df_SRR_complete_severe[
  with(df_SRR_complete_severe, order(-abundance)), ]

df_SRR_complete_severe$cumsum <- cumsum(
  df_SRR_complete_severe$abundance)
# table with abundant species
df_SRR_complete_severe <- subset(
  df_SRR_complete_severe, cumsum <= abund_spec)
df_SRR_complete_severe$cumsum <- NULL
df_SRR_complete_severe$abundance <- NULL
# extract absolute abundance values
df_SRR_severe_abundance_value <- rowMedians(
  as.matrix(df_SRR_complete_severe)) 
df_SRR_severe_abundance_species <- rownames(
  df_SRR_complete_severe)
df_SRR_severe_node_information <- data.frame(
  cbind(df_SRR_severe_abundance_value,
        df_SRR_severe_abundance_species))
df_SRR_severe_node_information <- separate(
  df_SRR_severe_node_information,
  col = 'df_SRR_severe_abundance_species',
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_SRR_severe_node_information$N <- NULL
# convert columns from factor to numeric
df_SRR_severe_node_information$df_SRR_severe_abundance_value <- as.numeric(
  as.character(df_SRR_severe_node_information$df_SRR_severe_abundance_value))
df_SRR_severe_node_information <- ddply(
  df_SRR_severe_node_information,"Genus",numcolwise(sum))
df_SRR_severe_node_information$df_SRR_severe_abundance_value <- round(
  df_SRR_severe_node_information$df_SRR_severe_abundance_value,2)
colnames(df_SRR_severe_node_information) <- c("Label", "Value")
df_SRR_complete_severe <- data.frame(
  t(df_SRR_complete_severe))

# generate correlation matrix (PI)
# severe
df_SRR_complete_severe_cor <- rcorr(
  as.matrix(df_SRR_complete_severe), 
  type = 'spearman')

# create node and edge lists
# severe
df_SRR_complete_severe_cor_Pval <- df_SRR_complete_severe_cor$P
df_SRR_complete_severe_COR <- df_SRR_complete_severe_cor$r

# generate edgelist
df_SRR_complete_severe_edges <- melt(
  df_SRR_complete_severe_cor_Pval)
df_SRR_complete_severe_coredges <- melt(
  df_SRR_complete_severe_COR)
df_SRR_complete_severe_edges$COR <- df_SRR_complete_severe_coredges$value
df_SRR_complete_severe_edges$value <- round(
  df_SRR_complete_severe_edges$value, 5)
df_SRR_complete_severe_edges$Label <- row.names(
  df_SRR_complete_severe_coredges)

colnames(df_SRR_complete_severe_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_SRR_complete_severe_edges) <- NULL
df_SRR_complete_severe_edges <- df_SRR_complete_severe_edges[
  complete.cases(df_SRR_complete_severe_edges), ]
df_SRR_complete_severe_edges_short <- subset(
  df_SRR_complete_severe_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_SRR_complete_severe_edges_short$cortype <- ifelse(
  df_SRR_complete_severe_edges_short$Weight > 0, "pos", "neg")
df_SRR_complete_severe_edges_short$Genus <- df_SRR_complete_severe_edges_short$Source
df_SRR_complete_severe_edges_short <- separate(
  df_SRR_complete_severe_edges_short,
  col = 'Genus', into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_SRR_complete_severe_edges_short$N <- NULL
rownames(df_SRR_complete_severe_edges_short) <- NULL
df_SRR_complete_severe_edges_short$Weight <- round(
  df_SRR_complete_severe_edges_short$Weight,2)
df_SRR_complete_severe_edges_short$Type <- "Undirected"

df_SRR_complete_severe_nodes <- select(
  df_SRR_complete_severe_edges_short, c(
    "Source", "Genus"))
df_SRR_complete_severe_nodes = df_SRR_complete_severe_nodes[
  !duplicated(df_SRR_complete_severe_nodes$Source),]
df_SRR_complete_severe_nodes <- data.frame(
  df_SRR_complete_severe_nodes)
rownames(df_SRR_complete_severe_nodes) <- NULL
colnames(df_SRR_complete_severe_nodes) <- c(
  "Id", "Label")
df_SRR_complete_severe_nodes1 <- merge(
  df_SRR_complete_severe_nodes, 
  df_SRR_severe_node_information[, c("Label", "Value")],
  by="Label")

# mild cases
# keep species that are non-zero in half of the patients
percent_SRR_complete_mild <- round(
  (ncol(df_SRR_complete_mild) * 33.3) / 100, 0)
df_SRR_complete_mild <- df_SRR_complete_mild[
  rowSums(df_SRR_complete_mild == 0) <= percent_SRR_complete_mild, ]

sum_SRR_mild = colSums(
  df_SRR_complete_mild)
sum_all_SRR_mild = sum(
  sum_SRR_mild)
df_SRR_complete_mild$abundance <- (rowSums(
  df_SRR_complete_mild[,1:ncol(df_SRR_complete_mild)]) / sum_all_SRR_mild) * 100
# sort abundance decreasing
df_SRR_complete_mild <- df_SRR_complete_mild[
  with(df_SRR_complete_mild,
       order(-abundance)), ]

df_SRR_complete_mild$cumsum <- cumsum(
  df_SRR_complete_mild$abundance)
# TABLE with abundant species
df_SRR_complete_mild <- subset(
  df_SRR_complete_mild, 
  cumsum <= abund_spec)
df_SRR_complete_mild$cumsum <- NULL
df_SRR_complete_mild$abundance <- NULL
# Extract absolute abundance values
df_SRR_mild_abundance_value <- rowMedians(
  as.matrix(df_SRR_complete_mild)) 
df_SRR_mild_abundance_species <- rownames(
  df_SRR_complete_mild)
df_SRR_mild_node_information <- data.frame(
  cbind(df_SRR_mild_abundance_value,
        df_SRR_mild_abundance_species))
df_SRR_mild_node_information <- separate(
  df_SRR_mild_node_information,
  col = 'df_SRR_mild_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_SRR_mild_node_information$N <- NULL
# convert columns from factor to numeric
df_SRR_mild_node_information$df_SRR_mild_abundance_value <- as.numeric(
  as.character(df_SRR_mild_node_information$df_SRR_mild_abundance_value))
df_SRR_mild_node_information <- ddply(
  df_SRR_mild_node_information,
  "Genus",
  numcolwise(sum))
df_SRR_mild_node_information$df_SRR_mild_abundance_value <- round(
  df_SRR_mild_node_information$df_SRR_mild_abundance_value,2)
colnames(df_SRR_mild_node_information) <- c("Label", "Value")
df_SRR_complete_mild <- data.frame(
  t(df_SRR_complete_mild))

# generate correlation matrix (PI)
# mild
df_SRR_complete_mild_cor <- rcorr(
  as.matrix(df_SRR_complete_mild), 
  type = 'spearman')
# create node and edge lists
# mild
df_SRR_complete_mild_cor_Pval <- df_SRR_complete_mild_cor$P
df_SRR_complete_mild_COR <- df_SRR_complete_mild_cor$r

# generate edgelist
df_SRR_complete_mild_edges <- melt(
  df_SRR_complete_mild_cor_Pval)
df_SRR_complete_mild_coredges <- melt(
  df_SRR_complete_mild_COR)
df_SRR_complete_mild_edges$COR <- df_SRR_complete_mild_coredges$value
df_SRR_complete_mild_edges$value <- round(
  df_SRR_complete_mild_edges$value, 5)
df_SRR_complete_mild_edges$Label <- row.names(
  df_SRR_complete_mild_coredges)

colnames(df_SRR_complete_mild_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_SRR_complete_mild_edges) <- NULL
df_SRR_complete_mild_edges <- df_SRR_complete_mild_edges[
  complete.cases(df_SRR_complete_mild_edges), ]
df_SRR_complete_mild_edges_short <- subset(
  df_SRR_complete_mild_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_SRR_complete_mild_edges_short$cortype <- ifelse(
  df_SRR_complete_mild_edges_short$Weight > 0, "pos", "neg")
df_SRR_complete_mild_edges_short$Genus <- df_SRR_complete_mild_edges_short$Source
df_SRR_complete_mild_edges_short <- separate(
  df_SRR_complete_mild_edges_short, 
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_SRR_complete_mild_edges_short$N <- NULL
rownames(df_SRR_complete_mild_edges_short) <- NULL
df_SRR_complete_mild_edges_short$Weight <- round(
  df_SRR_complete_mild_edges_short$Weight,2)
df_SRR_complete_mild_edges_short$Type <- "Undirected"

df_SRR_complete_mild_nodes <- select(
  df_SRR_complete_mild_edges_short, 
  c("Source", "Genus"))
df_SRR_complete_mild_nodes = df_SRR_complete_mild_nodes[
  !duplicated(df_SRR_complete_mild_nodes$Source),]
df_SRR_complete_mild_nodes <- data.frame(
  df_SRR_complete_mild_nodes)
rownames(df_SRR_complete_mild_nodes) <- NULL
colnames(df_SRR_complete_mild_nodes) <- c("Id", "Label")
df_SRR_complete_mild_nodes1 <- merge(
  df_SRR_complete_mild_nodes,
  df_SRR_mild_node_information[, c("Label", "Value")], by="Label")

######################################################################
# solid data ####
# import data file
df_solid_complete <- read_delim(
  "input_files/solid_CF_mixed.csv", 
  ";", escape_double = FALSE, trim_ws = TRUE)
df_solid_complete <- data.frame(
  df_solid_complete)
df_solid_complete <- df_solid_complete[
  complete.cases(df_solid_complete), ]
df_solid_complete$Species <- str_replace(
  df_solid_complete$Species, " ", "_")
df_solid_complete <- ddply(
  df_solid_complete,"Species",numcolwise(sum))

# remove viral/phage/fungi entries
df_solid_complete <- df_solid_complete[
  - grep("phage", df_solid_complete$Species),]
df_solid_complete <- df_solid_complete[
  - grep("virus", df_solid_complete$Species),]
df_solid_complete <- df_solid_complete[
  - grep("albicans", df_solid_complete$Species),]
df_solid_complete <- df_solid_complete[
  - grep("fumigatus", df_solid_complete$Species),]
df_solid_complete <- df_solid_complete[
  - grep("dubliniensis", df_solid_complete$Species),]
df_solid_complete <- df_solid_complete[
  - grep("pinophilus", df_solid_complete$Species),]
rownames(df_solid_complete) <- df_solid_complete$Species
df_solid_complete$Species <- NULL
# separate PS/PI entries
df_solid_PS <- df_solid_complete[
  , grep("PS", colnames(df_solid_complete))]
df_solid_PI <- df_solid_complete[
  , grep("PI", colnames(df_solid_complete))]
# separate by sample type
# nasal lavage
df_solid_PS_N <- df_solid_PS[
  , grep("_N.", colnames(df_solid_PS))]
n_PS_N = ncol(df_solid_PS_N)
df_solid_PI_N <- df_solid_PI[
  , grep("_N.", colnames(df_solid_PI))]
n_PI_N = ncol(df_solid_PI_N)
# sputum
df_solid_PS_S <- df_solid_PS[
  , grep("_S.", colnames(df_solid_PS))]
n_PS_S = ncol(df_solid_PS_S)
df_solid_PI_S <- df_solid_PI[
  , grep("_S.", colnames(df_solid_PI))]
n_PI_S = ncol(df_solid_PI_S)
# throat
df_solid_PS_T <- df_solid_PS[
  , grep("_T.", colnames(df_solid_PS))]
n_PS_T = ncol(df_solid_PS_T)
df_solid_PI_T <- df_solid_PI[
  , grep("_T.", colnames(df_solid_PI))]
n_PI_T = ncol(df_solid_PI_T)

# separate by age group 
# nasal lavage
# A
df_solid_PS_N_A <- df_solid_PS_N[
  , grep("PSA", colnames(df_solid_PS_N))]
n_PS_N_A = ncol(df_solid_PS_N_A)
df_solid_PI_N_A <- df_solid_PI_N[
  , grep("PIA", colnames(df_solid_PI_N))]
n_PI_N_A = ncol(df_solid_PI_N_A)

# B
df_solid_PS_N_B <- df_solid_PS_N[
  , grep("PSB", colnames(df_solid_PS_N))]
n_PS_N_B = ncol(df_solid_PS_N_B)
df_solid_PI_N_B <- df_solid_PI_N[
  , grep("PIB", colnames(df_solid_PI_N))]
n_PI_N_B = ncol(df_solid_PI_N_B)

# C
df_solid_PS_N_C <- df_solid_PS_N[
  , grep("PSC", colnames(df_solid_PS_N))]
n_PS_N_C = ncol(df_solid_PS_N_C)

df_solid_PI_N_C <- df_solid_PI_N[
  , grep("PIC", colnames(df_solid_PI_N))]
n_PI_N_C = ncol(df_solid_PI_N_C)

# sputum
# A
df_solid_PS_S_A <- df_solid_PS_S[
  , grep("PSA", colnames(df_solid_PS_S))]
n_PS_S_A = ncol(df_solid_PS_S_A)
df_solid_PI_S_A <- df_solid_PI_S[
  , grep("PIA", colnames(df_solid_PI_S))]
n_PI_S_A = ncol(df_solid_PI_S_A)

# B
df_solid_PS_S_B <- df_solid_PS_S[
  , grep("PSA", colnames(df_solid_PS_S))]
n_PS_S_B = ncol(df_solid_PS_S_B)

df_solid_PI_S_B <- df_solid_PI_S[
  , grep("PIB", colnames(df_solid_PI_S))]
n_PI_S_B = ncol(df_solid_PI_S_B)

# C
df_solid_PS_S_C <- df_solid_PS_S[
  , grep("PSC", colnames(df_solid_PS_S))]
n_PS_S_C = ncol(df_solid_PS_S_C)

df_solid_PI_S_C <- df_solid_PI_S[
  , grep("PIC", colnames(df_solid_PI_S))]
n_PI_S_C = ncol(df_solid_PI_S_C)

# throat
# A
df_solid_PS_T_A <- df_solid_PS_T[
  , grep("PSA", colnames(df_solid_PS_T))]
n_PS_T_A = ncol(df_solid_PS_T_A)

df_solid_PI_T_A <- df_solid_PI_T[
  , grep("PIA", colnames(df_solid_PI_T))]
n_PI_T_A = ncol(df_solid_PI_T_A)

# B
df_solid_PS_T_B <- df_solid_PS_T[
  , grep("PSB", colnames(df_solid_PS_T))]
n_PS_T_B = ncol(df_solid_PS_T_B)

df_solid_PI_T_B <- df_solid_PI_T[
  , grep("PIB", colnames(df_solid_PI_T))]
n_PI_T_B = ncol(df_solid_PI_T_B)

# C
df_solid_PS_T_C <- df_solid_PS_T[
  , grep("PSC", colnames(df_solid_PS_T))]
n_PS_T_C = ncol(df_solid_PS_T_C)

df_solid_PI_T_C <- df_solid_PI_T[
  , grep("PIC", colnames(df_solid_PI_T))]
n_PI_T_C = ncol(df_solid_PI_T_C)

sample = c(
  "sputum", "sputum_A", "sputum_B", "sputum_C",
  "nasal", "nasal_A", "nasal_B", "nasal_C",
  "throat", "throat_A", "throat_B", "throat_C")
col_pi = c(
  n_PI_S, n_PI_S_A,n_PI_S_B, n_PI_S_C,
  n_PI_N, n_PI_N_A,n_PI_N_B, n_PI_N_C,
  n_PI_T, n_PI_T_A,n_PI_T_B, n_PI_T_C)
col_ps = c(
  n_PS_S, n_PS_S_A,n_PS_S_B, n_PS_S_C,
  n_PS_N, n_PS_N_A,n_PS_N_B, n_PS_N_C,
           n_PS_T, n_PS_T_A,n_PS_T_B, n_PS_T_C)
df_n <-data.frame(
  cbind(sample, col_pi, col_ps))

# continue with sputum
# PS, A
sum_solid_PS_S_A = colSums(
  df_solid_PS_S_A)
sum_all_solid_PS_S_A = sum(
  sum_solid_PS_S_A)
df_solid_PS_S_A$abundance <- (rowSums(
  df_solid_PS_S_A[,1:ncol(df_solid_PS_S_A)]) / sum_all_solid_PS_S_A) * 100
# sort abundance decreasing
df_solid_PS_S_A <- df_solid_PS_S_A[
  with(df_solid_PS_S_A, order(-abundance)), ]

df_solid_PS_S_A$cumsum <- cumsum(
  df_solid_PS_S_A$abundance)
# table with abundant species
df_solid_PS_S_A <- subset(
  df_solid_PS_S_A, cumsum <= abund_spec)
df_solid_PS_S_A$cumsum <- NULL
df_solid_PS_S_A$abundance <- NULL


# extract absolute abundance values
df_solid_PS_S_A_abundance_value <- rowMedians(
  as.matrix(df_solid_PS_S_A)) 
df_solid_PS_S_A_abundance_species <- rownames(
  df_solid_PS_S_A)
df_solid_PS_S_A_node_information <- data.frame(
  cbind(df_solid_PS_S_A_abundance_value,
        df_solid_PS_S_A_abundance_species))
df_solid_PS_S_A_node_information <- separate(
  df_solid_PS_S_A_node_information,
  col = 'df_solid_PS_S_A_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_solid_PS_S_A_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PS_S_A_node_information$df_solid_PS_S_A_abundance_value <- as.numeric(
  as.character(df_solid_PS_S_A_node_information$df_solid_PS_S_A_abundance_value))
df_solid_PS_S_A_node_information <- ddply(
  df_solid_PS_S_A_node_information,
  "Genus",
  numcolwise(sum))
df_solid_PS_S_A_node_information$df_solid_PS_S_A_abundance_value <- round(
  df_solid_PS_S_A_node_information$df_solid_PS_S_A_abundance_value,2)
colnames(df_solid_PS_S_A_node_information) <- c("Label", "Value")
df_solid_PS_S_A <- data.frame(
  t(df_solid_PS_S_A))

# generate correlation matrix (PI)
# PS_S_A
df_solid_PS_S_A_cor <- rcorr(
  as.matrix(df_solid_PS_S_A), 
  type = 'spearman')

# create node and edge lists
# PS_S_A
df_solid_PS_S_A_cor_Pval <- df_solid_PS_S_A_cor$P
ddf_solid_PS_S_A_COR <- df_solid_PS_S_A_cor$r

# generate edgelist
df_solid_PS_S_A_edges <- melt(
  df_solid_PS_S_A_cor_Pval)
df_solid_PS_S_A_coredges <- melt(
  ddf_solid_PS_S_A_COR)
df_solid_PS_S_A_edges$COR <- df_solid_PS_S_A_coredges$value
df_solid_PS_S_A_edges$value <- round(
  df_solid_PS_S_A_edges$value, 5)
df_solid_PS_S_A_edges$Label <- row.names(
  df_solid_PS_S_A_coredges)

colnames(df_solid_PS_S_A_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PS_S_A_edges) <- NULL
df_solid_PS_S_A_edges <- df_solid_PS_S_A_edges[
  complete.cases(df_solid_PS_S_A_edges), ]
df_solid_PS_S_A_edges_short <- subset(
  df_solid_PS_S_A_edges, 
  pValue < sig_level & Weight >= weight_val |
    pValue < sig_level & Weight <= weight_val2)
df_solid_PS_S_A_edges_short$cortype <-ifelse(
  df_solid_PS_S_A_edges_short$Weight > 0, "pos", "neg")
df_solid_PS_S_A_edges_short$Genus <- df_solid_PS_S_A_edges_short$Source
df_solid_PS_S_A_edges_short <- separate(
  df_solid_PS_S_A_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PS_S_A_edges_short$N <- NULL
rownames(df_solid_PS_S_A_edges_short) <- NULL
df_solid_PS_S_A_edges_short$Weight <- round(
  df_solid_PS_S_A_edges_short$Weight,2)
df_solid_PS_S_A_edges_short$Type <- "Undirected"

df_solid_PS_S_A_edges_nodes <- select(
  df_solid_PS_S_A_edges_short, 
  c("Source", "Genus"))
df_solid_PS_S_A_edges_nodes = df_solid_PS_S_A_edges_nodes[
  !duplicated(df_solid_PS_S_A_edges_nodes$Source),]
df_solid_PS_S_A_edges_nodes <- data.frame(
  df_solid_PS_S_A_edges_nodes)
rownames(df_solid_PS_S_A_edges_nodes) <- NULL
colnames(df_solid_PS_S_A_edges_nodes) <- c("Id", "Label")
df_solid_PS_S_A_edges_nodes1 <- merge(
  df_solid_PS_S_A_edges_nodes, 
  df_solid_PS_S_A_node_information[
    , c("Label", "Value")], 
  by="Label")

# Continue with sputum
# PS, B
sum_solid_PS_S_B = colSums(
  df_solid_PS_S_B)
sum_all_solid_PS_S_B = sum(
  sum_solid_PS_S_B)
df_solid_PS_S_B$abundance <- (
  rowSums(df_solid_PS_S_B[,1:ncol(df_solid_PS_S_B)]) / sum_all_solid_PS_S_B) * 100
# sort abundance decreasing
df_solid_PS_S_B <- df_solid_PS_S_B[
  with(df_solid_PS_S_B, order(-abundance)), ]

df_solid_PS_S_B$cumsum <- cumsum(
  df_solid_PS_S_B$abundance)
# TABLE with abundant species
df_solid_PS_S_B <- subset(
  df_solid_PS_S_B, cumsum <= abund_spec)
df_solid_PS_S_B$cumsum <- NULL
df_solid_PS_S_B$abundance <- NULL

# Extract absolute abundance values
df_solid_PS_S_B_abundance_value <- rowMedians(
  as.matrix(df_solid_PS_S_B)) 
df_solid_PS_S_B_abundance_species <- rownames(df_solid_PS_S_B)
df_solid_PS_S_B_node_information <- data.frame(
  cbind(df_solid_PS_S_B_abundance_value,
        df_solid_PS_S_B_abundance_species))
df_solid_PS_S_B_node_information <- separate(
  df_solid_PS_S_B_node_information,
  col = 'df_solid_PS_S_B_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PS_S_B_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PS_S_B_node_information$df_solid_PS_S_B_abundance_value <- as.numeric(
  as.character(df_solid_PS_S_B_node_information$df_solid_PS_S_B_abundance_value))
df_solid_PS_S_B_node_information <- ddply(
  df_solid_PS_S_B_node_information, 
  "Genus",
  numcolwise(sum))
df_solid_PS_S_B_node_information$df_solid_PS_S_B_abundance_value <- round(
  df_solid_PS_S_B_node_information$df_solid_PS_S_B_abundance_value,2)
colnames(df_solid_PS_S_B_node_information) <- c("Label", "Value")
df_solid_PS_S_B <- data.frame(
  t(df_solid_PS_S_B))

# Generate correlation matrix (PI)
# PS_S_B
df_solid_PS_S_B_cor <- rcorr(
  as.matrix(df_solid_PS_S_B), 
  type = 'spearman')

# Create node and edge lists
# PS_S_B
df_solid_PS_S_B_cor_Pval <- df_solid_PS_S_B_cor$P
ddf_solid_PS_S_B_COR <- df_solid_PS_S_B_cor$r

# generate edgelist
df_solid_PS_S_B_edges <- melt(
  df_solid_PS_S_B_cor_Pval)
df_solid_PS_S_B_coredges <- melt(
  ddf_solid_PS_S_B_COR)
df_solid_PS_S_B_edges$COR <- df_solid_PS_S_B_coredges$value
df_solid_PS_S_B_edges$value <- round(
  df_solid_PS_S_B_edges$value, 5)
df_solid_PS_S_B_edges$Label <- row.names(df_solid_PS_S_B_coredges)

colnames(df_solid_PS_S_B_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PS_S_B_edges) <- NULL
df_solid_PS_S_B_edges <- df_solid_PS_S_B_edges[
  complete.cases(df_solid_PS_S_B_edges), ]
df_solid_PS_S_B_edges_short <- subset(
  df_solid_PS_S_B_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PS_S_B_edges_short$cortype <- ifelse(
  df_solid_PS_S_B_edges_short$Weight > 0, "pos", "neg")
df_solid_PS_S_B_edges_short$Genus <- df_solid_PS_S_B_edges_short$Source
df_solid_PS_S_B_edges_short <- separate(
  df_solid_PS_S_B_edges_short,
  col = 'Genus', into=c('Genus', 'N'),
  sep='_', remove=TRUE)
df_solid_PS_S_B_edges_short$N <- NULL
rownames(df_solid_PS_S_B_edges_short) <- NULL
df_solid_PS_S_B_edges_short$Weight <- round(
  df_solid_PS_S_B_edges_short$Weight,2)
df_solid_PS_S_B_edges_short$Type <- "Undirected"

df_solid_PS_S_B_edges_nodes <- select(
  df_solid_PS_S_B_edges_short, 
  c("Source", "Genus"))
df_solid_PS_S_B_edges_nodes = df_solid_PS_S_B_edges_nodes[
  !duplicated(df_solid_PS_S_B_edges_nodes$Source),]
df_solid_PS_S_B_edges_nodes <- data.frame(
  df_solid_PS_S_B_edges_nodes)
rownames(df_solid_PS_S_B_edges_nodes) <- NULL
colnames(df_solid_PS_S_B_edges_nodes) <- c("Id", "Label")
df_solid_PS_S_B_edges_nodes1 <- merge(
  df_solid_PS_S_B_edges_nodes, 
  df_solid_PS_S_B_node_information[
    , c("Label", "Value")], 
  by="Label")


# Continue with sputum
# PS, C
sum_solid_PS_S_C = colSums(
  df_solid_PS_S_C)
sum_all_solid_PS_S_C = sum(
  sum_solid_PS_S_C)
df_solid_PS_S_C$abundance <- (rowSums(
  df_solid_PS_S_C[,1:ncol(df_solid_PS_S_C)]) / sum_all_solid_PS_S_C) * 100
# sort abundance decreasing
df_solid_PS_S_C <- df_solid_PS_S_C[
  with(df_solid_PS_S_C, order(-abundance)), ]

df_solid_PS_S_C$cumsum <- cumsum(
  df_solid_PS_S_C$abundance)
# table with abundant species
df_solid_PS_S_C <- subset(
  df_solid_PS_S_C, cumsum <= abund_spec)
df_solid_PS_S_C$cumsum <- NULL
df_solid_PS_S_C$abundance <- NULL

# Extract absolute abundance values
df_solid_PS_S_C_abundance_value <- rowMedians(
  as.matrix(df_solid_PS_S_C)) 
df_solid_PS_S_C_abundance_species <- rownames(
  df_solid_PS_S_C)
df_solid_PS_S_C_node_information <- data.frame(
  cbind(df_solid_PS_S_C_abundance_value,
        df_solid_PS_S_C_abundance_species))
df_solid_PS_S_C_node_information <- separate(
  df_solid_PS_S_C_node_information,
  col = 'df_solid_PS_S_C_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PS_S_C_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PS_S_C_node_information$df_solid_PS_S_C_abundance_value <- as.numeric(
  as.character(df_solid_PS_S_C_node_information$df_solid_PS_S_C_abundance_value))
df_solid_PS_S_C_node_information <- ddply(
  df_solid_PS_S_C_node_information,
  "Genus",
  numcolwise(sum))
df_solid_PS_S_C_node_information$df_solid_PS_S_C_abundance_value <- round(
  df_solid_PS_S_C_node_information$df_solid_PS_S_C_abundance_value,2)
colnames(df_solid_PS_S_C_node_information) <- c("Label", "Value")
df_solid_PS_S_C <- data.frame(
  t(df_solid_PS_S_C))

# Generate correlation matrix (PI)
# PS_S_C
df_solid_PS_S_C_cor <- rcorr(
  as.matrix(df_solid_PS_S_C), 
  type = 'spearman')

# Create node and edge lists
# PS_S_C
df_solid_PS_S_C_cor_Pval <- df_solid_PS_S_C_cor$P
ddf_solid_PS_S_C_COR <- df_solid_PS_S_C_cor$r

# generate edgelist
df_solid_PS_S_C_edges <- melt(
  df_solid_PS_S_C_cor_Pval)
df_solid_PS_S_C_coredges <- melt(
  ddf_solid_PS_S_C_COR)
df_solid_PS_S_C_edges$COR <- df_solid_PS_S_C_coredges$value
df_solid_PS_S_C_edges$value <- round(
  df_solid_PS_S_C_edges$value, 5)
df_solid_PS_S_C_edges$Label <- row.names(df_solid_PS_S_C_coredges)

colnames(df_solid_PS_S_C_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PS_S_C_edges) <- NULL
df_solid_PS_S_C_edges <- df_solid_PS_S_C_edges[
  complete.cases(df_solid_PS_S_C_edges), ]
df_solid_PS_S_C_edges_short <- subset(
  df_solid_PS_S_C_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PS_S_C_edges_short$cortype <- ifelse(
  df_solid_PS_S_C_edges_short$Weight > 0, "pos", "neg")
df_solid_PS_S_C_edges_short$Genus <- df_solid_PS_S_C_edges_short$Source
df_solid_PS_S_C_edges_short <- separate(
  df_solid_PS_S_C_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PS_S_C_edges_short$N <- NULL
rownames(df_solid_PS_S_C_edges_short) <- NULL
df_solid_PS_S_C_edges_short$Weight <- round(
  df_solid_PS_S_C_edges_short$Weight,2)
df_solid_PS_S_C_edges_short$Type <- "Undirected"

df_solid_PS_S_C_edges_nodes <- select(
  df_solid_PS_S_C_edges_short, 
  c("Source", "Genus"))
df_solid_PS_S_C_edges_nodes = df_solid_PS_S_C_edges_nodes[
  !duplicated(df_solid_PS_S_C_edges_nodes$Source),]
df_solid_PS_S_C_edges_nodes <- data.frame(
  df_solid_PS_S_C_edges_nodes)
rownames(df_solid_PS_S_C_edges_nodes) <- NULL
colnames(df_solid_PS_S_C_edges_nodes) <- c("Id", "Label")
df_solid_PS_S_C_edges_nodes1 <- merge(
  df_solid_PS_S_C_edges_nodes, 
  df_solid_PS_S_C_node_information[
    , c("Label", "Value")], by="Label")

# Continue with sputum
# PI, A
sum_solid_PI_S_A = colSums(
  df_solid_PI_S_A)
sum_all_solid_PI_S_A = sum(
  sum_solid_PI_S_A)
df_solid_PI_S_A$abundance <- (rowSums(
  df_solid_PI_S_A[,1:ncol(df_solid_PI_S_A)]) / sum_all_solid_PI_S_A) * 100
# sort abundance decreasing
df_solid_PI_S_A <- df_solid_PI_S_A[
  with(df_solid_PI_S_A, order(-abundance)), ]

df_solid_PI_S_A$cumsum <- cumsum(
  df_solid_PI_S_A$abundance)
# table with abundant species
df_solid_PI_S_A <- subset(
  df_solid_PI_S_A, cumsum <= abund_spec)
df_solid_PI_S_A$cumsum <- NULL
df_solid_PI_S_A$abundance <- NULL

# Extract absolute abundance values
df_solid_PI_S_A_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_S_A)) 
df_solid_PI_S_A_abundance_species <- rownames(df_solid_PI_S_A)
df_solid_PI_S_A_node_information <- data.frame(
  cbind(df_solid_PI_S_A_abundance_value,
        df_solid_PI_S_A_abundance_species))
df_solid_PI_S_A_node_information <- separate(
  df_solid_PI_S_A_node_information,
  col = 'df_solid_PI_S_A_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_S_A_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_S_A_node_information$df_solid_PI_S_A_abundance_value <- as.numeric(
  as.character(df_solid_PI_S_A_node_information$df_solid_PI_S_A_abundance_value))
df_solid_PI_S_A_node_information <- ddply(
  df_solid_PI_S_A_node_information,"Genus",numcolwise(sum))
df_solid_PI_S_A_node_information$df_solid_PI_S_A_abundance_value <- round(
  df_solid_PI_S_A_node_information$df_solid_PI_S_A_abundance_value,2)
colnames(df_solid_PI_S_A_node_information) <- c("Label", "Value")
df_solid_PI_S_A <- data.frame(
  t(df_solid_PI_S_A))

# Generate correlation matrix (PI)
# PI_S_A
df_solid_PI_S_A_cor <- rcorr(
  as.matrix(df_solid_PI_S_A), 
  type = 'spearman')

# Create node and edge lists
# PI_S_A
df_solid_PI_S_A_cor_Pval <- df_solid_PI_S_A_cor$P
ddf_solid_PI_S_A_COR <- df_solid_PI_S_A_cor$r

# generate edgelist
df_solid_PI_S_A_edges <- melt(
  df_solid_PI_S_A_cor_Pval)
df_solid_PI_S_A_coredges <- melt(
  ddf_solid_PI_S_A_COR)
df_solid_PI_S_A_edges$COR <- df_solid_PI_S_A_coredges$value
df_solid_PI_S_A_edges$value <- round(
  df_solid_PI_S_A_edges$value, 5)
df_solid_PI_S_A_edges$Label <- row.names(
  df_solid_PI_S_A_coredges)

colnames(df_solid_PI_S_A_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_S_A_edges) <- NULL
df_solid_PI_S_A_edges <- df_solid_PI_S_A_edges[
  complete.cases(df_solid_PI_S_A_edges), ]
df_solid_PI_S_A_edges_short <- subset(
  df_solid_PI_S_A_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_S_A_edges_short$cortype <- ifelse(
  df_solid_PI_S_A_edges_short$Weight > 0, "pos", "neg")
df_solid_PI_S_A_edges_short$Genus <- df_solid_PI_S_A_edges_short$Source
df_solid_PI_S_A_edges_short <- separate(
  df_solid_PI_S_A_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_S_A_edges_short$N <- NULL
rownames(df_solid_PI_S_A_edges_short) <- NULL
df_solid_PI_S_A_edges_short$Weight <- round(
  df_solid_PI_S_A_edges_short$Weight,2)
df_solid_PI_S_A_edges_short$Type <- "Undirected"

df_solid_PI_S_A_edges_nodes <- select(
  df_solid_PI_S_A_edges_short, 
  c("Source", "Genus"))
df_solid_PI_S_A_edges_nodes = df_solid_PI_S_A_edges_nodes[
  !duplicated(df_solid_PI_S_A_edges_nodes$Source),]
df_solid_PI_S_A_edges_nodes <- data.frame(
  df_solid_PI_S_A_edges_nodes)
rownames(df_solid_PI_S_A_edges_nodes) <- NULL
colnames(df_solid_PI_S_A_edges_nodes) <- c("Id", "Label")
df_solid_PI_S_A_edges_nodes1 <- merge(
  df_solid_PI_S_A_edges_nodes, 
  df_solid_PI_S_A_node_information[
    , c("Label", "Value")], 
  by="Label")

# PI, B
sum_solid_PI_S_B = colSums(
  df_solid_PI_S_B)
sum_all_solid_PI_S_B = sum(
  sum_solid_PI_S_B)
df_solid_PI_S_B$abundance <- (rowSums(
  df_solid_PI_S_B[,1:ncol(df_solid_PI_S_B)]) / sum_all_solid_PI_S_B) * 100
# sort abundance decreasing
df_solid_PI_S_B <- df_solid_PI_S_B[
  with(df_solid_PI_S_B, order(-abundance)), ]

df_solid_PI_S_B$cumsum <- cumsum(
  df_solid_PI_S_B$abundance)
# table with abundant species
df_solid_PI_S_B <- subset(
  df_solid_PI_S_B, cumsum <= abund_spec)
df_solid_PI_S_B$cumsum <- NULL
df_solid_PI_S_B$abundance <- NULL

# extract absolute abundance values
df_solid_PI_S_B_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_S_B)) 
df_solid_PI_S_B_abundance_species <- rownames(df_solid_PI_S_B)
df_solid_PI_S_B_node_information <- data.frame(
  cbind(df_solid_PI_S_B_abundance_value,
        df_solid_PI_S_B_abundance_species))
df_solid_PI_S_B_node_information <- separate(
  df_solid_PI_S_B_node_information,
  col = 'df_solid_PI_S_B_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_S_B_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_S_B_node_information$df_solid_PI_S_B_abundance_value <- as.numeric(
  as.character(df_solid_PI_S_B_node_information$df_solid_PI_S_B_abundance_value))
df_solid_PI_S_B_node_information <- ddply(
  df_solid_PI_S_B_node_information,"Genus",numcolwise(sum))
df_solid_PI_S_B_node_information$df_solid_PI_S_B_abundance_value <- round(
  df_solid_PI_S_B_node_information$df_solid_PI_S_B_abundance_value,2)
colnames(df_solid_PI_S_B_node_information) <- c("Label", "Value")
df_solid_PI_S_B <- data.frame(
  t(df_solid_PI_S_B))

# generate correlation matrix (PI)
# PI_S_B
df_solid_PI_S_B_cor <- rcorr(
  as.matrix(df_solid_PI_S_B), 
  type = 'spearman')

# create node and edge lists
# PI_S_B
df_solid_PI_S_B_cor_Pval <- df_solid_PI_S_B_cor$P
ddf_solid_PI_S_B_COR <- df_solid_PI_S_B_cor$r

# generate edgelist
df_solid_PI_S_B_edges <- melt(
  df_solid_PI_S_B_cor_Pval)
df_solid_PI_S_B_coredges <- melt(
  ddf_solid_PI_S_B_COR)
df_solid_PI_S_B_edges$COR <- df_solid_PI_S_B_coredges$value
df_solid_PI_S_B_edges$value <- round(
  df_solid_PI_S_B_edges$value, 5)
df_solid_PI_S_B_edges$Label <- row.names(
  df_solid_PI_S_B_coredges)

colnames(df_solid_PI_S_B_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_S_B_edges) <- NULL
df_solid_PI_S_B_edges <- df_solid_PI_S_B_edges[
  complete.cases(df_solid_PI_S_B_edges), ]
df_solid_PI_S_B_edges_short <- subset(
  df_solid_PI_S_B_edges, pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_S_B_edges_short$cortype <- ifelse(
  df_solid_PI_S_B_edges_short$Weight > 0, "pos", "neg")

df_solid_PI_S_B_edges_short$Genus <- df_solid_PI_S_B_edges_short$Source
df_solid_PI_S_B_edges_short <- separate(
  df_solid_PI_S_B_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'),
  sep='_', 
  remove=TRUE)
df_solid_PI_S_B_edges_short$N <- NULL
rownames(df_solid_PI_S_B_edges_short) <- NULL
df_solid_PI_S_B_edges_short$Weight <- round(
  df_solid_PI_S_B_edges_short$Weight,2)
df_solid_PI_S_B_edges_short$Type <- "Undirected"

df_solid_PI_S_B_edges_nodes <- select(
  df_solid_PI_S_B_edges_short, 
  c("Source", "Genus"))
df_solid_PI_S_B_edges_nodes = df_solid_PI_S_B_edges_nodes[
  !duplicated(df_solid_PI_S_B_edges_nodes$Source),]
df_solid_PI_S_B_edges_nodes <- data.frame(
  df_solid_PI_S_B_edges_nodes)
rownames(df_solid_PI_S_B_edges_nodes) <- NULL
colnames(df_solid_PI_S_B_edges_nodes) <- c("Id", "Label")
df_solid_PI_S_B_edges_nodes1 <-merge(
  df_solid_PI_S_B_edges_nodes, 
  df_solid_PI_S_B_node_information[
    , c("Label", "Value")], 
  by="Label")


# PI, C
sum_solid_PI_S_C = colSums(
  df_solid_PI_S_C)
sum_all_solid_PI_S_C = sum(
  sum_solid_PI_S_C)
df_solid_PI_S_C$abundance <- (rowSums(
  df_solid_PI_S_C[,1:ncol(df_solid_PI_S_C)]) / sum_all_solid_PI_S_C) * 100
# sort abundance decreasing
df_solid_PI_S_C <- df_solid_PI_S_C[
  with(df_solid_PI_S_C, order(-abundance)), ]
df_solid_PI_S_C$cumsum <- cumsum(
  df_solid_PI_S_C$abundance)
# table with abundant species
df_solid_PI_S_C <- subset(
  df_solid_PI_S_C, cumsum <= abund_spec)
df_solid_PI_S_C$cumsum <- NULL
df_solid_PI_S_C$abundance <- NULL

# extract absolute abundance values
df_solid_PI_S_C_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_S_C)) 
df_solid_PI_S_C_abundance_species <- rownames(
  df_solid_PI_S_C)
df_solid_PI_S_C_node_information <- data.frame(
  cbind(df_solid_PI_S_C_abundance_value,
        df_solid_PI_S_C_abundance_species))
df_solid_PI_S_C_node_information <- separate(
  df_solid_PI_S_C_node_information,
  col = 'df_solid_PI_S_C_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_S_C_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_S_C_node_information$df_solid_PI_S_C_abundance_value <- as.numeric(
  as.character(df_solid_PI_S_C_node_information$df_solid_PI_S_C_abundance_value))
df_solid_PI_S_C_node_information <- ddply(
  df_solid_PI_S_C_node_information,
  "Genus",
  numcolwise(sum))
df_solid_PI_S_C_node_information$df_solid_PI_S_C_abundance_value <- round(
  df_solid_PI_S_C_node_information$df_solid_PI_S_C_abundance_value,2)
colnames(df_solid_PI_S_C_node_information) <- c("Label", "Value")
df_solid_PI_S_C <- data.frame(
  t(df_solid_PI_S_C))

# generate correlation matrix (PI)
# PI_S_C
df_solid_PI_S_C_cor <- rcorr(
  as.matrix(df_solid_PI_S_C), 
  type = 'spearman')

# create node and edge lists
# PI_S_C
df_solid_PI_S_C_cor_Pval <- df_solid_PI_S_C_cor$P
ddf_solid_PI_S_C_COR <- df_solid_PI_S_C_cor$r

# generate edgelist
df_solid_PI_S_C_edges <- melt(
  df_solid_PI_S_C_cor_Pval)
df_solid_PI_S_C_coredges <- melt(
  ddf_solid_PI_S_C_COR)
df_solid_PI_S_C_edges$COR <- df_solid_PI_S_C_coredges$value
df_solid_PI_S_C_edges$value <- round(
  df_solid_PI_S_C_edges$value, 5)
df_solid_PI_S_C_edges$Label <- row.names(
  df_solid_PI_S_C_coredges)

colnames(df_solid_PI_S_C_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_S_C_edges) <- NULL
df_solid_PI_S_C_edges <- df_solid_PI_S_C_edges[
  complete.cases(df_solid_PI_S_C_edges), ]
df_solid_PI_S_C_edges_short <- subset(
  df_solid_PI_S_C_edges, 
  pValue < sig_level & Weight >= weight_val |
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_S_C_edges_short$cortype <- ifelse(
  df_solid_PI_S_C_edges_short$Weight > 0, "pos", "neg")
df_solid_PI_S_C_edges_short$Genus <- df_solid_PI_S_C_edges_short$Source
df_solid_PI_S_C_edges_short <- separate(
  df_solid_PI_S_C_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_solid_PI_S_C_edges_short$N <- NULL
rownames(df_solid_PI_S_C_edges_short) <- NULL
df_solid_PI_S_C_edges_short$Weight <- round(
  df_solid_PI_S_C_edges_short$Weight,2)
df_solid_PI_S_C_edges_short$Type <- "Undirected"

df_solid_PI_S_C_edges_nodes <- select(
  df_solid_PI_S_C_edges_short, c("Source", "Genus"))
df_solid_PI_S_C_edges_nodes = df_solid_PI_S_C_edges_nodes[
  !duplicated(df_solid_PI_S_C_edges_nodes$Source),]
df_solid_PI_S_C_edges_nodes <- data.frame(
  df_solid_PI_S_C_edges_nodes)
rownames(df_solid_PI_S_C_edges_nodes) <- NULL
colnames(df_solid_PI_S_C_edges_nodes) <- c("Id", "Label")
df_solid_PI_S_C_edges_nodes1 <- merge(
  df_solid_PI_S_C_edges_nodes, 
  df_solid_PI_S_C_node_information[
    ,c("Label", "Value")], 
  by="Label")


# continue with nasal lavage
# PS, A
sum_solid_PS_N_A = colSums(
  df_solid_PS_N_A)
sum_all_solid_PS_N_A = sum(
  sum_solid_PS_N_A)
df_solid_PS_N_A$abundance <- (rowSums(
  df_solid_PS_N_A[,1:ncol(df_solid_PS_N_A)]) / sum_all_solid_PS_N_A) * 100
# sort abundance decreasing
df_solid_PS_N_A <- df_solid_PS_N_A[
  with(df_solid_PS_N_A, order(-abundance)), ]

df_solid_PS_N_A$cumsum <- cumsum(
  df_solid_PS_N_A$abundance)
# table with abundant species
df_solid_PS_N_A <- subset(
  df_solid_PS_N_A, cumsum <= abund_spec)
df_solid_PS_N_A$cumsum <- NULL
df_solid_PS_N_A$abundance <- NULL


# extract absolute abundance values
df_solid_PS_N_A_abundance_value <- rowMedians(
  as.matrix(df_solid_PS_N_A)) 
df_solid_PS_N_A_abundance_species <- rownames(
  df_solid_PS_N_A)
df_solid_PS_N_A_node_information <- data.frame(
  cbind(df_solid_PS_N_A_abundance_value,
        df_solid_PS_N_A_abundance_species))
df_solid_PS_N_A_node_information <- separate(
  df_solid_PS_N_A_node_information,
  col = 'df_solid_PS_N_A_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_solid_PS_N_A_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PS_N_A_node_information$df_solid_PS_N_A_abundance_value <- as.numeric(
  as.character(df_solid_PS_N_A_node_information$df_solid_PS_N_A_abundance_value))
df_solid_PS_N_A_node_information <- ddply(
  df_solid_PS_N_A_node_information,
  "Genus",
  numcolwise(sum))
df_solid_PS_N_A_node_information$df_solid_PS_N_A_abundance_value <- round(
  df_solid_PS_N_A_node_information$df_solid_PS_N_A_abundance_value,2)
colnames(df_solid_PS_N_A_node_information) <- c("Label", "Value")
df_solid_PS_N_A <- data.frame(
  t(df_solid_PS_N_A))

# generate correlation matrix (PI)
# PS_N_A
df_solid_PS_N_A_cor <- rcorr(
  as.matrix(df_solid_PS_N_A), 
  type = 'spearman')

# create node and edge lists
# PS_N_A
df_solid_PS_N_A_cor_Pval <- df_solid_PS_N_A_cor$P
ddf_solid_PS_N_A_COR <- df_solid_PS_N_A_cor$r

# generate edgelist
df_solid_PS_N_A_edges <- melt(
  df_solid_PS_N_A_cor_Pval)
df_solid_PS_N_A_coredges <- melt(
  ddf_solid_PS_N_A_COR)
df_solid_PS_N_A_edges$COR <- df_solid_PS_N_A_coredges$value
df_solid_PS_N_A_edges$value <- round(
  df_solid_PS_N_A_edges$value, 5)
df_solid_PS_N_A_edges$Label <- row.names(
  df_solid_PS_N_A_coredges)

colnames(df_solid_PS_N_A_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PS_N_A_edges) <- NULL
df_solid_PS_N_A_edges <- df_solid_PS_N_A_edges[
  complete.cases(df_solid_PS_N_A_edges), ]
df_solid_PS_N_A_edges_short <- subset(
  df_solid_PS_N_A_edges, 
  pValue < sig_level & Weight >= weight_val |
    pValue < sig_level & Weight <= weight_val2)
df_solid_PS_N_A_edges_short$cortype <-ifelse(
  df_solid_PS_N_A_edges_short$Weight > 0, "pos", "neg")
df_solid_PS_N_A_edges_short$Genus <- df_solid_PS_N_A_edges_short$Source
df_solid_PS_N_A_edges_short <- separate(
  df_solid_PS_N_A_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PS_N_A_edges_short$N <- NULL
rownames(df_solid_PS_N_A_edges_short) <- NULL
df_solid_PS_N_A_edges_short$Weight <- round(
  df_solid_PS_N_A_edges_short$Weight,2)
df_solid_PS_N_A_edges_short$Type <- "Undirected"

df_solid_PS_N_A_edges_nodes <- select(
  df_solid_PS_N_A_edges_short, 
  c("Source", "Genus"))
df_solid_PS_N_A_edges_nodes = df_solid_PS_N_A_edges_nodes[
  !duplicated(df_solid_PS_N_A_edges_nodes$Source),]
df_solid_PS_N_A_edges_nodes <- data.frame(
  df_solid_PS_N_A_edges_nodes)
rownames(df_solid_PS_N_A_edges_nodes) <- NULL
colnames(df_solid_PS_N_A_edges_nodes) <- c("Id", "Label")
df_solid_PS_N_A_edges_nodes1 <- merge(
  df_solid_PS_N_A_edges_nodes, 
  df_solid_PS_N_A_node_information[
    , c("Label", "Value")], 
  by="Label")

# Continue with sputum
# PS, B (skip, less than 4 variables)
# PS, C (skip, less than 4 variables)
# PI, A
sum_solid_PI_N_A = colSums(
  df_solid_PI_N_A)
sum_all_solid_PI_N_A = sum(
  sum_solid_PI_N_A)
df_solid_PI_N_A$abundance <- (rowSums(
  df_solid_PI_N_A[,1:ncol(df_solid_PI_N_A)]) / sum_all_solid_PI_N_A) * 100
# sort abundance decreasing
df_solid_PI_N_A <- df_solid_PI_N_A[
  with(df_solid_PI_N_A, order(-abundance)), ]

df_solid_PI_N_A$cumsum <- cumsum(
  df_solid_PI_N_A$abundance)
# table with abundant species
df_solid_PI_N_A <- subset(
  df_solid_PI_N_A, cumsum <= abund_spec)
df_solid_PI_N_A$cumsum <- NULL
df_solid_PI_N_A$abundance <- NULL

# Extract absolute abundance values
df_solid_PI_N_A_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_N_A)) 
df_solid_PI_N_A_abundance_species <- rownames(df_solid_PI_N_A)
df_solid_PI_N_A_node_information <- data.frame(
  cbind(df_solid_PI_N_A_abundance_value,
        df_solid_PI_N_A_abundance_species))
df_solid_PI_N_A_node_information <- separate(
  df_solid_PI_N_A_node_information,
  col = 'df_solid_PI_N_A_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_N_A_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_N_A_node_information$df_solid_PI_N_A_abundance_value <- as.numeric(
  as.character(df_solid_PI_N_A_node_information$df_solid_PI_N_A_abundance_value))
df_solid_PI_N_A_node_information <- ddply(
  df_solid_PI_N_A_node_information,"Genus",numcolwise(sum))
df_solid_PI_N_A_node_information$df_solid_PI_N_A_abundance_value <- round(
  df_solid_PI_N_A_node_information$df_solid_PI_N_A_abundance_value,2)
colnames(df_solid_PI_N_A_node_information) <- c("Label", "Value")
df_solid_PI_N_A <- data.frame(
  t(df_solid_PI_N_A))

# Generate correlation matrix (PI)
# PI_N_A
df_solid_PI_N_A_cor <- rcorr(
  as.matrix(df_solid_PI_N_A), 
  type = 'spearman')

# Create node and edge lists
# PI_N_A
df_solid_PI_N_A_cor_Pval <- df_solid_PI_N_A_cor$P
ddf_solid_PI_N_A_COR <- df_solid_PI_N_A_cor$r

# generate edgelist
df_solid_PI_N_A_edges <- melt(
  df_solid_PI_N_A_cor_Pval)
df_solid_PI_N_A_coredges <- melt(
  ddf_solid_PI_N_A_COR)
df_solid_PI_N_A_edges$COR <- df_solid_PI_N_A_coredges$value
df_solid_PI_N_A_edges$value <- round(
  df_solid_PI_N_A_edges$value, 5)
df_solid_PI_N_A_edges$Label <- row.names(
  df_solid_PI_N_A_coredges)

colnames(df_solid_PI_N_A_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_N_A_edges) <- NULL
df_solid_PI_N_A_edges <- df_solid_PI_N_A_edges[
  complete.cases(df_solid_PI_N_A_edges), ]
df_solid_PI_N_A_edges_short <- subset(
  df_solid_PI_N_A_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_N_A_edges_short$cortype <- ifelse(
  df_solid_PI_N_A_edges_short$Weight > 0, "pos", "neg")
df_solid_PI_N_A_edges_short$Genus <- df_solid_PI_N_A_edges_short$Source
df_solid_PI_N_A_edges_short <- separate(
  df_solid_PI_N_A_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_N_A_edges_short$N <- NULL
rownames(df_solid_PI_N_A_edges_short) <- NULL
df_solid_PI_N_A_edges_short$Weight <- round(
  df_solid_PI_N_A_edges_short$Weight,2)
df_solid_PI_N_A_edges_short$Type <- "Undirected"

df_solid_PI_N_A_edges_nodes <- select(
  df_solid_PI_N_A_edges_short, 
  c("Source", "Genus"))
df_solid_PI_N_A_edges_nodes = df_solid_PI_N_A_edges_nodes[
  !duplicated(df_solid_PI_N_A_edges_nodes$Source),]
df_solid_PI_N_A_edges_nodes <- data.frame(
  df_solid_PI_N_A_edges_nodes)
rownames(df_solid_PI_N_A_edges_nodes) <- NULL
colnames(df_solid_PI_N_A_edges_nodes) <- c("Id", "Label")
df_solid_PI_N_A_edges_nodes1 <- merge(
  df_solid_PI_N_A_edges_nodes, 
  df_solid_PI_N_A_node_information[
    , c("Label", "Value")], 
  by="Label")

# PI, B
sum_solid_PI_N_B = colSums(
  df_solid_PI_N_B)
sum_all_solid_PI_N_B = sum(
  sum_solid_PI_N_B)
df_solid_PI_N_B$abundance <- (rowSums(
  df_solid_PI_N_B[,1:ncol(df_solid_PI_N_B)]) / sum_all_solid_PI_N_B) * 100
# sort abundance decreasing
df_solid_PI_N_B <- df_solid_PI_N_B[
  with(df_solid_PI_N_B, order(-abundance)), ]

df_solid_PI_N_B$cumsum <- cumsum(
  df_solid_PI_N_B$abundance)
# table with abundant species
df_solid_PI_N_B <- subset(
  df_solid_PI_N_B, cumsum <= abund_spec)
df_solid_PI_N_B$cumsum <- NULL
df_solid_PI_N_B$abundance <- NULL

# extract absolute abundance values
df_solid_PI_N_B_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_N_B)) 
df_solid_PI_N_B_abundance_species <- rownames(df_solid_PI_N_B)
df_solid_PI_N_B_node_information <- data.frame(
  cbind(df_solid_PI_N_B_abundance_value,
        df_solid_PI_N_B_abundance_species))
df_solid_PI_N_B_node_information <- separate(
  df_solid_PI_N_B_node_information,
  col = 'df_solid_PI_N_B_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_N_B_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_N_B_node_information$df_solid_PI_N_B_abundance_value <- as.numeric(
  as.character(df_solid_PI_N_B_node_information$df_solid_PI_N_B_abundance_value))
df_solid_PI_N_B_node_information <- ddply(
  df_solid_PI_N_B_node_information,"Genus",numcolwise(sum))
df_solid_PI_N_B_node_information$df_solid_PI_N_B_abundance_value <- round(
  df_solid_PI_N_B_node_information$df_solid_PI_N_B_abundance_value,2)
colnames(df_solid_PI_N_B_node_information) <- c("Label", "Value")
df_solid_PI_N_B <- data.frame(
  t(df_solid_PI_N_B))

# generate correlation matrix (PI)
# PI_N_B
df_solid_PI_N_B_cor <- rcorr(
  as.matrix(df_solid_PI_N_B), 
  type = 'spearman')

# create node and edge lists
# PI_N_B
df_solid_PI_N_B_cor_Pval <- df_solid_PI_N_B_cor$P
ddf_solid_PI_N_B_COR <- df_solid_PI_N_B_cor$r

# generate edgelist
df_solid_PI_N_B_edges <- melt(
  df_solid_PI_N_B_cor_Pval)
df_solid_PI_N_B_coredges <- melt(
  ddf_solid_PI_N_B_COR)
df_solid_PI_N_B_edges$COR <- df_solid_PI_N_B_coredges$value
df_solid_PI_N_B_edges$value <- round(
  df_solid_PI_N_B_edges$value, 5)
df_solid_PI_N_B_edges$Label <- row.names(
  df_solid_PI_N_B_coredges)

colnames(df_solid_PI_N_B_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_N_B_edges) <- NULL
df_solid_PI_N_B_edges <- df_solid_PI_N_B_edges[
  complete.cases(df_solid_PI_N_B_edges), ]
df_solid_PI_N_B_edges_short <- subset(
  df_solid_PI_N_B_edges, pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_N_B_edges_short$cortype <- ifelse(
  df_solid_PI_N_B_edges_short$Weight > 0, "pos", "neg")

df_solid_PI_N_B_edges_short$Genus <- df_solid_PI_N_B_edges_short$Source
df_solid_PI_N_B_edges_short <- separate(
  df_solid_PI_N_B_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'),
  sep='_', 
  remove=TRUE)
df_solid_PI_N_B_edges_short$N <- NULL
rownames(df_solid_PI_N_B_edges_short) <- NULL
df_solid_PI_N_B_edges_short$Weight <- round(
  df_solid_PI_N_B_edges_short$Weight,2)
df_solid_PI_N_B_edges_short$Type <- "Undirected"

df_solid_PI_N_B_edges_nodes <- select(
  df_solid_PI_N_B_edges_short, 
  c("Source", "Genus"))
df_solid_PI_N_B_edges_nodes = df_solid_PI_N_B_edges_nodes[
  !duplicated(df_solid_PI_N_B_edges_nodes$Source),]
df_solid_PI_N_B_edges_nodes <- data.frame(
  df_solid_PI_N_B_edges_nodes)
rownames(df_solid_PI_N_B_edges_nodes) <- NULL
colnames(df_solid_PI_N_B_edges_nodes) <- c("Id", "Label")
df_solid_PI_N_B_edges_nodes1 <-merge(
  df_solid_PI_N_B_edges_nodes, 
  df_solid_PI_N_B_node_information[
    , c("Label", "Value")], 
  by="Label")


# PI, C
sum_solid_PI_N_C = colSums(
  df_solid_PI_N_C)
sum_all_solid_PI_N_C = sum(
  sum_solid_PI_N_C)
df_solid_PI_N_C$abundance <- (rowSums(
  df_solid_PI_N_C[,1:ncol(df_solid_PI_N_C)]) / sum_all_solid_PI_N_C) * 100
# sort abundance decreasing
df_solid_PI_N_C <- df_solid_PI_N_C[
  with(df_solid_PI_N_C, order(-abundance)), ]
df_solid_PI_N_C$cumsum <- cumsum(
  df_solid_PI_N_C$abundance)
# table with abundant species
df_solid_PI_N_C <- subset(
  df_solid_PI_N_C, cumsum <= abund_spec)
df_solid_PI_N_C$cumsum <- NULL
df_solid_PI_N_C$abundance <- NULL

# extract absolute abundance values
df_solid_PI_N_C_abundance_value <- rowMedians(
  as.matrix(df_solid_PI_N_C)) 
df_solid_PI_N_C_abundance_species <- rownames(
  df_solid_PI_N_C)
df_solid_PI_N_C_node_information <- data.frame(
  cbind(df_solid_PI_N_C_abundance_value,
        df_solid_PI_N_C_abundance_species))
df_solid_PI_N_C_node_information <- separate(
  df_solid_PI_N_C_node_information,
  col = 'df_solid_PI_N_C_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_solid_PI_N_C_node_information$N <- NULL
# convert columns from factor to numeric
df_solid_PI_N_C_node_information$df_solid_PI_N_C_abundance_value <- as.numeric(
  as.character(df_solid_PI_N_C_node_information$df_solid_PI_N_C_abundance_value))
df_solid_PI_N_C_node_information <- ddply(
  df_solid_PI_N_C_node_information,
  "Genus",
  numcolwise(sum))
df_solid_PI_N_C_node_information$df_solid_PI_N_C_abundance_value <- round(
  df_solid_PI_N_C_node_information$df_solid_PI_N_C_abundance_value,2)
colnames(df_solid_PI_N_C_node_information) <- c("Label", "Value")
df_solid_PI_N_C <- data.frame(
  t(df_solid_PI_N_C))

# generate correlation matrix (PI)
# PI_N_C
df_solid_PI_N_C_cor <- rcorr(
  as.matrix(df_solid_PI_N_C), 
  type = 'spearman')

# create node and edge lists
# PI_N_C
df_solid_PI_N_C_cor_Pval <- df_solid_PI_N_C_cor$P
ddf_solid_PI_N_C_COR <- df_solid_PI_N_C_cor$r

# generate edgelist
df_solid_PI_N_C_edges <- melt(
  df_solid_PI_N_C_cor_Pval)
df_solid_PI_N_C_coredges <- melt(
  ddf_solid_PI_N_C_COR)
df_solid_PI_N_C_edges$COR <- df_solid_PI_N_C_coredges$value
df_solid_PI_N_C_edges$value <- round(
  df_solid_PI_N_C_edges$value, 5)
df_solid_PI_N_C_edges$Label <- row.names(
  df_solid_PI_N_C_coredges)

colnames(df_solid_PI_N_C_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_solid_PI_N_C_edges) <- NULL
df_solid_PI_N_C_edges <- df_solid_PI_N_C_edges[
  complete.cases(df_solid_PI_N_C_edges), ]
df_solid_PI_N_C_edges_short <- subset(
  df_solid_PI_N_C_edges, 
  pValue < sig_level & Weight >= weight_val |
    pValue < sig_level & Weight <= weight_val2)
df_solid_PI_N_C_edges_short$cortype <- ifelse(
  df_solid_PI_N_C_edges_short$Weight > 0, "pos", "neg")
df_solid_PI_N_C_edges_short$Genus <- df_solid_PI_N_C_edges_short$Source
df_solid_PI_N_C_edges_short <- separate(
  df_solid_PI_N_C_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_solid_PI_N_C_edges_short$N <- NULL
rownames(df_solid_PI_N_C_edges_short) <- NULL
df_solid_PI_N_C_edges_short$Weight <- round(
  df_solid_PI_N_C_edges_short$Weight,2)
df_solid_PI_N_C_edges_short$Type <- "Undirected"

df_solid_PI_N_C_edges_nodes <- select(
  df_solid_PI_N_C_edges_short, c("Source", "Genus"))
df_solid_PI_N_C_edges_nodes = df_solid_PI_N_C_edges_nodes[
  !duplicated(df_solid_PI_N_C_edges_nodes$Source),]
df_solid_PI_N_C_edges_nodes <- data.frame(
  df_solid_PI_N_C_edges_nodes)
rownames(df_solid_PI_N_C_edges_nodes) <- NULL
colnames(df_solid_PI_N_C_edges_nodes) <- c("Id", "Label")
df_solid_PI_N_C_edges_nodes1 <- merge(
  df_solid_PI_N_C_edges_nodes, 
  df_solid_PI_N_C_node_information[
    ,c("Label", "Value")], 
  by="Label")


########################################################################
# Newborns ####
df_newborns_complete <- read_delim(
  "input_files/illumina_newborns_coughSwabs.csv", ";", 
  escape_double = FALSE, trim_ws = TRUE)
df_newborns_complete <- data.frame(
  df_newborns_complete)
df_newborns_complete <- df_newborns_complete[
  complete.cases(df_newborns_complete), ]
df_newborns_complete$Species <- str_replace(
  df_newborns_complete$Species, " ", "_")
df_newborns_complete <- ddply(
  df_newborns_complete,"Species",numcolwise(sum))
df_newborns_complete <- df_newborns_complete[
  - grep("phage", df_newborns_complete$Species),]
df_newborns_complete <- df_newborns_complete[
  - grep("virus", df_newborns_complete$Species),]
df_newborns_complete <- df_newborns_complete[
  - grep("_sp", df_newborns_complete$Species),]
rownames(df_newborns_complete) <- df_newborns_complete$Species
df_newborns_complete$Species <- NULL

# select CF and healthy pre school childen
healthy_pre <- c(
  "KGCF01", "KGCF05", "KGCF06", "KGCF19", "KGCF21", 
  "KGCF25", "KGCF27", "KGCF28", "KGCF29", "KGCF30",
  "KGCF33", "KGCF35", "KGCF38", "KGCF40")

cf_pre <- c(
  "EMCF01","EMCF03","EMCF04","EMCF17","EMCF19",
  "EMCF21","EMCF22","EMCF25","EMCF26","EMCF30",
  "EMCF32","EMCF33","EMCF34","EMCF35","EMCF36",
  "EMCF37")

df_healthy_pre <- select(
  df_newborns_complete, 
  healthy_pre)

df_CF_pre <- select(
  df_newborns_complete, 
  cf_pre)

# continue with healthy preschoolers
sum_healthy_pre = colSums(
  df_healthy_pre)
sum_all_healthy_pre = sum(
  sum_healthy_pre)
df_healthy_pre$abundance <- (rowSums(
  df_healthy_pre[,1:ncol(df_healthy_pre)]) / sum_all_healthy_pre) * 100
# sort abundance decreasing
df_healthy_pre <- df_healthy_pre[
    with(df_healthy_pre, order(-abundance)), ]

df_healthy_pre$cumsum <- cumsum(
  df_healthy_pre$abundance)
# table with abundant species
df_healthy_pre <- subset(
  df_healthy_pre, 
  cumsum <= abund_spec)
df_healthy_pre$cumsum <- NULL
df_healthy_pre$abundance <- NULL

# extract absolute abundance values
df_healthy_pre_abundance_value <- rowMedians(
  as.matrix(df_healthy_pre)) 
df_healthy_pre_abundance_species <- rownames(
  df_healthy_pre)
df_healthy_pre_node_information <- data.frame(
  cbind(df_healthy_pre_abundance_value,
        df_healthy_pre_abundance_species))
df_healthy_pre_node_information <- separate(
  df_healthy_pre_node_information,
  col = 'df_healthy_pre_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_healthy_pre_node_information$N <- NULL
# convert columns from factor to numeric
df_healthy_pre_node_information$df_healthy_pre_abundance_value <- as.numeric(
  as.character(df_healthy_pre_node_information$df_healthy_pre_abundance_value))
df_healthy_pre_node_information <- ddply(
  df_healthy_pre_node_information,
  "Genus",
  numcolwise(sum))
df_healthy_pre_node_information$df_healthy_pre_abundance_value <- round(
  df_healthy_pre_node_information$df_healthy_pre_abundance_value,2)
colnames(df_healthy_pre_node_information) <- c("Label", "Value")
df_healthy_pre <- data.frame(
  t(df_healthy_pre))

# generate correlation matrix 
df_healthy_pre_cor <- rcorr(
  as.matrix(df_healthy_pre), 
  type = 'spearman')

# create node and edge lists
df_healthy_pre_cor_Pval <- df_healthy_pre_cor$P
df_healthy_pre_COR <- df_healthy_pre_cor$r

# generate edgelist
df_healthy_pre_edges <- melt(
  df_healthy_pre_cor_Pval)
df_healthy_pre_coredges <- melt(
  df_healthy_pre_COR)
df_healthy_pre_edges$COR <- df_healthy_pre_coredges$value
df_healthy_pre_edges$value <- round(
  df_healthy_pre_edges$value, 5)
df_healthy_pre_edges$Label <- row.names(
  df_healthy_pre_coredges)

colnames(df_healthy_pre_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_healthy_pre_edges) <- NULL
df_healthy_pre_edges <- df_healthy_pre_edges[
  complete.cases(df_healthy_pre_edges), ]
df_healthy_pre_edges_short <- subset(
  df_healthy_pre_edges, pValue < sig_level & Weight >= weight_val |
    pValue < sig_level & Weight <= weight_val2)
df_healthy_pre_edges_short$cortype <- ifelse(
  df_healthy_pre_edges_short$Weight > 0, "pos", "neg")
df_healthy_pre_edges_short$Genus <- df_healthy_pre_edges_short$Source
df_healthy_pre_edges_short <- separate(
  df_healthy_pre_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_healthy_pre_edges_short$N <- NULL
rownames(df_healthy_pre_edges_short) <- NULL
df_healthy_pre_edges_short$Weight <- round(
  df_healthy_pre_edges_short$Weight,2)
df_healthy_pre_edges_short$Type <- "Undirected"

df_healthy_pre_edges_nodes <- select(
  df_healthy_pre_edges_short, 
  c("Source", "Genus"))
df_healthy_pre_edges_nodes = df_healthy_pre_edges_nodes[
  !duplicated(df_healthy_pre_edges_nodes$Source),]
df_healthy_pre_edges_nodes <- data.frame(
  df_healthy_pre_edges_nodes)
rownames(df_healthy_pre_edges_nodes) <- NULL
colnames(df_healthy_pre_edges_nodes) <- c("Id", "Label")
df_healthy_pre_edges_nodes1 <- merge(
  df_healthy_pre_edges_nodes, 
  df_healthy_pre_node_information[, c("Label", "Value")], by="Label")

# continue with CF preschoolers
sum_CF_pre = colSums(
  df_CF_pre)
sum_all_CF_pre = sum(
  sum_CF_pre)
df_CF_pre$abundance <- (rowSums(
  df_CF_pre[,1:ncol(df_CF_pre)]) / sum_all_CF_pre) * 100
# sort abundance decreasing
df_CF_pre <- df_CF_pre[
  with(df_CF_pre, order(-abundance)), ]

df_CF_pre$cumsum <- cumsum(
  df_CF_pre$abundance)
# table with abundant species
df_CF_pre <- subset(
  df_CF_pre, cumsum <= abund_spec)
df_CF_pre$cumsum <- NULL
df_CF_pre$abundance <- NULL

# extract absolute abundance values
df_CF_pre_abundance_value <- rowMedians(
  as.matrix(df_CF_pre)) 
df_CF_pre_abundance_species <- rownames(df_CF_pre)
df_CF_pre_node_information <- data.frame(
  cbind(df_CF_pre_abundance_value,
        df_CF_pre_abundance_species))
df_CF_pre_node_information <- separate(
  df_CF_pre_node_information,
  col = 'df_CF_pre_abundance_species', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_CF_pre_node_information$N <- NULL
# convert columns from factor to numeric
df_CF_pre_node_information$df_CF_pre_abundance_value <- as.numeric(
  as.character(df_CF_pre_node_information$df_CF_pre_abundance_value))
df_CF_pre_node_information <- ddply(
  df_CF_pre_node_information,
  "Genus",
  numcolwise(sum))
df_CF_pre_node_information$df_CF_pre_abundance_value <- round(
  df_CF_pre_node_information$df_CF_pre_abundance_value,2)
colnames(df_CF_pre_node_information) <- c("Label", "Value")
df_CF_pre <- data.frame(
  t(df_CF_pre))

# generate correlation matrix (PI)
df_CF_pre_cor <- rcorr(
  as.matrix(df_CF_pre), 
  type = 'spearman')

# create node and edge lists
df_CF_pre_cor_Pval <- df_CF_pre_cor$P
df_CF_pre_COR <- df_CF_pre_cor$r

# generate edgelist
df_CF_pre_edges <- melt(
  df_CF_pre_cor_Pval)
df_CF_pre_coredges <- melt(
  df_CF_pre_COR)
df_CF_pre_edges$COR <- df_CF_pre_coredges$value
df_CF_pre_edges$value <- round(
  df_CF_pre_edges$value, 5)
df_CF_pre_edges$Label <- row.names(
  df_CF_pre_coredges)

colnames(df_CF_pre_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_CF_pre_edges) <- NULL
df_CF_pre_edges <- df_CF_pre_edges[
  complete.cases(df_CF_pre_edges), ]
df_CF_pre_edges_short <- subset(
  df_CF_pre_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_CF_pre_edges_short$cortype <-ifelse(
  df_CF_pre_edges_short$Weight > 0, "pos", "neg")
df_CF_pre_edges_short$Genus <- df_CF_pre_edges_short$Source
df_CF_pre_edges_short <- separate(
  df_CF_pre_edges_short,
  col = 'Genus', 
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_CF_pre_edges_short$N <- NULL
rownames(df_CF_pre_edges_short) <- NULL
df_CF_pre_edges_short$Weight <- round(
  df_CF_pre_edges_short$Weight,2)
df_CF_pre_edges_short$Type <- "Undirected"

df_CF_pre_edges_nodes <- select(
  df_CF_pre_edges_short, 
  c("Source", "Genus"))
df_CF_pre_edges_nodes = df_CF_pre_edges_nodes[
  !duplicated(df_CF_pre_edges_nodes$Source),]
df_CF_pre_edges_nodes <- data.frame(
  df_CF_pre_edges_nodes)
rownames(df_CF_pre_edges_nodes) <- NULL
colnames(df_CF_pre_edges_nodes) <- c("Id", "Label")
df_CF_pre_edges_nodes1 <- merge(
  df_CF_pre_edges_nodes, 
  df_CF_pre_node_information[
    , c("Label", "Value")], 
  by="Label")

# Healthy sputum ####
#Smoker and non-smokers 
healthy_sputum <- read_delim(
  "input_files/smoker_nonsmoker.csv",
  ";", escape_double = FALSE, trim_ws = TRUE)
healthy_sputum <- data.frame(
  healthy_sputum)

healthy_sputum <- ddply(
  healthy_sputum,"Species", numcolwise(sum))
healthy_sputum <- healthy_sputum[
  - grep("phage", healthy_sputum$Species),]
healthy_sputum <- healthy_sputum[
  - grep("virus", healthy_sputum$Species),]
healthy_sputum <- healthy_sputum[
  - grep("albicans", healthy_sputum$Species),]

healthy_sputum <- healthy_sputum[
  complete.cases(healthy_sputum), ]
healthy_sputum$Species <- str_replace(
  healthy_sputum$Species, " ", "_")
rownames(healthy_sputum) <- healthy_sputum$Species
healthy_sputum$Species <- NULL

colnames(healthy_sputum) <- c(
  "10_SPU_SMO",	"11_SPU_SMO",
  "12_SPU_SMO", "13_SPU_NO_SMO",
  "14_SPU_SMO", "15_SPU_NO_SMO",	
  "17_SPU_SMO", "18_SPU_SMO",
  "19_SPU_NO_SMO", "20_SPU_SMO",
  "22_SPU_NO_SMO", "2_SPU_SMO",
  "3_SPU_NO_SMO", "4_SPU_SMO", 
  "9_SPU_SMO")

healthy_sputum_smo <- c(
  "10_SPU_SMO", "11_SPU_SMO", "12_SPU_SMO", 
  "14_SPU_SMO", "17_SPU_SMO", "18_SPU_SMO", "20_SPU_SMO", 
  "2_SPU_SMO", "4_SPU_SMO", "9_SPU_SMO")
healthy_sputum_nonsmo <- c(
  "13_SPU_NO_SMO", "15_SPU_NO_SMO", "19_SPU_NO_SMO", 
  "22_SPU_NO_SMO", "3_SPU_NO_SMO")

df_healthy_sputum_smo <- select(
  healthy_sputum, healthy_sputum_smo)

df_healthy_sputum_no_smo <- select(
  healthy_sputum, healthy_sputum_nonsmo)

# smokers
sum_healthy_sputum_smo = colSums(
  df_healthy_sputum_smo)
sum_all_healthy_sputum_smo = sum(
  sum_healthy_sputum_smo)
df_healthy_sputum_smo$abundance <- (rowSums(
  df_healthy_sputum_smo[,1:ncol(df_healthy_sputum_smo)]) / sum_all_healthy_sputum_smo) * 100
# sort abundance decreasing
df_healthy_sputum_smo <- df_healthy_sputum_smo[
  with(df_healthy_sputum_smo, order(-abundance)), ]

df_healthy_sputum_smo$cumsum <- cumsum(
  df_healthy_sputum_smo$abundance)
# table with abundant species
df_healthy_sputum_smo <- subset(
  df_healthy_sputum_smo, cumsum <= abund_spec)
df_healthy_sputum_smo$cumsum <- NULL
df_healthy_sputum_smo$abundance <- NULL
# extract absolute abundance values
df_healthy_sputum_smo_abundance_value <- rowMedians(
  as.matrix(df_healthy_sputum_smo)) 
df_healthy_sputum_smo_abundance_species <- rownames(
  df_healthy_sputum_smo)
df_healthy_sputum_smo_information <- data.frame(
  cbind(df_healthy_sputum_smo_abundance_value,
        df_healthy_sputum_smo_abundance_species))
df_healthy_sputum_smo_information <- separate(
  df_healthy_sputum_smo_information,
  col = 'df_healthy_sputum_smo_abundance_species',
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_healthy_sputum_smo_information$N <- NULL
# convert columns from factor to numeric
df_healthy_sputum_smo_information$df_healthy_sputum_smo_abundance_value <- as.numeric(
  as.character(df_healthy_sputum_smo_information$df_healthy_sputum_smo_abundance_value))
df_healthy_sputum_smo_information <- ddply(
  df_healthy_sputum_smo_information,"Genus",numcolwise(sum))
df_healthy_sputum_smo_information$df_healthy_sputum_smo_abundance_value <- round(
  df_healthy_sputum_smo_information$df_healthy_sputum_smo_abundance_value,2)
colnames(df_healthy_sputum_smo_information) <- c("Label", "Value")
df_healthy_sputum_smo <- data.frame(
  t(df_healthy_sputum_smo))

# generate correlation matrix (PI)
# severe
df_healthy_sputum_smo_cor <- rcorr(
  as.matrix(df_healthy_sputum_smo), 
  type = 'spearman')

# create node and edge lists
# severe
df_healthy_sputum_smo_cor_Pval <- df_healthy_sputum_smo_cor$P
df_healthy_sputum_smo_COR <- df_healthy_sputum_smo_cor$r

# generate edgelist
df_healthy_sputum_smo_edges <- melt(
  df_healthy_sputum_smo_cor_Pval)
df_healthy_sputum_smo_coredges <- melt(
  df_healthy_sputum_smo_COR)
df_healthy_sputum_smo_edges$COR <- df_healthy_sputum_smo_coredges$value
df_healthy_sputum_smo_edges$value <- round(
  df_healthy_sputum_smo_edges$value, 5)
df_healthy_sputum_smo_edges$Label <- row.names(
  df_healthy_sputum_smo_coredges)

colnames(df_healthy_sputum_smo_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_healthy_sputum_smo_edges) <- NULL
df_healthy_sputum_smo_edges <- df_healthy_sputum_smo_edges[
  complete.cases(df_healthy_sputum_smo_edges), ]
df_healthy_sputum_smo_edges_short <- subset(
  df_healthy_sputum_smo_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_healthy_sputum_smo_edges_short$cortype <- ifelse(
  df_healthy_sputum_smo_edges_short$Weight > 0, "pos", "neg")
df_healthy_sputum_smo_edges_short$Genus <- df_healthy_sputum_smo_edges_short$Source
df_healthy_sputum_smo_edges_short <- separate(
  df_healthy_sputum_smo_edges_short,
  col = 'Genus', into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_healthy_sputum_smo_edges_short$N <- NULL
rownames(df_healthy_sputum_smo_edges_short) <- NULL
df_healthy_sputum_smo_edges_short$Weight <- round(
  df_healthy_sputum_smo_edges_short$Weight,2)
df_healthy_sputum_smo_edges_short$Type <- "Undirected"

df_healthy_sputum_smo_nodes <- select(
  df_healthy_sputum_smo_edges_short, c(
    "Source", "Genus"))
df_healthy_sputum_smo_nodes = df_healthy_sputum_smo_nodes[
  !duplicated(df_healthy_sputum_smo_nodes$Source),]
df_healthy_sputum_smo_nodes <- data.frame(
  df_healthy_sputum_smo_nodes)
rownames(df_healthy_sputum_smo_nodes) <- NULL
colnames(df_healthy_sputum_smo_nodes) <- c(
  "Id", "Label")
df_healthy_sputum_smo_nodes1 <- merge(
  df_healthy_sputum_smo_nodes, 
  df_healthy_sputum_smo_information[, c("Label", "Value")],
  by="Label")

# Non-smokers ###
df_healthy_sputum_no_smo <- select(
  healthy_sputum, healthy_sputum_nonsmo)

# non smokers
sum_healthy_sputum_no_smo = colSums(
  df_healthy_sputum_no_smo)
sum_all_healthy_sputum_no_smo = sum(
  sum_healthy_sputum_no_smo)
df_healthy_sputum_no_smo$abundance <- (rowSums(
  df_healthy_sputum_no_smo[,1:ncol(df_healthy_sputum_no_smo)]) / sum_all_healthy_sputum_no_smo) * 100
# sort abundance decreasing
df_healthy_sputum_no_smo <- df_healthy_sputum_no_smo[
  with(df_healthy_sputum_no_smo, order(-abundance)), ]

df_healthy_sputum_no_smo$cumsum <- cumsum(
  df_healthy_sputum_no_smo$abundance)
# table with abundant species
df_healthy_sputum_no_smo <- subset(
  df_healthy_sputum_no_smo, cumsum <= abund_spec)
df_healthy_sputum_no_smo$cumsum <- NULL
df_healthy_sputum_no_smo$abundance <- NULL
# extract absolute abundance values
df_healthy_sputum_no_smo_abundance_value <- rowMedians(
  as.matrix(df_healthy_sputum_no_smo)) 
df_healthy_sputum_no_smo_abundance_species <- rownames(
  df_healthy_sputum_no_smo)
df_healthy_sputum_no_smo_information <- data.frame(
  cbind(df_healthy_sputum_no_smo_abundance_value,
        df_healthy_sputum_no_smo_abundance_species))
df_healthy_sputum_no_smo_information <- separate(
  df_healthy_sputum_no_smo_information,
  col = 'df_healthy_sputum_no_smo_abundance_species',
  into=c('Genus', 'N'), 
  sep='_', 
  remove=TRUE)
df_healthy_sputum_no_smo_information$N <- NULL
# convert columns from factor to numeric
df_healthy_sputum_no_smo_information$df_healthy_sputum_no_smo_abundance_value <- as.numeric(
  as.character(df_healthy_sputum_no_smo_information$df_healthy_sputum_no_smo_abundance_value))
df_healthy_sputum_no_smo_information <- ddply(
  df_healthy_sputum_no_smo_information,"Genus",numcolwise(sum))
df_healthy_sputum_no_smo_information$df_healthy_sputum_no_smo_abundance_value <- round(
  df_healthy_sputum_no_smo_information$df_healthy_sputum_no_smo_abundance_value,2)
colnames(df_healthy_sputum_no_smo_information) <- c("Label", "Value")
df_healthy_sputum_no_smo <- data.frame(
  t(df_healthy_sputum_no_smo))

# generate correlation matrix (PI)
# severe
df_healthy_sputum_no_smo_cor <- rcorr(
  as.matrix(df_healthy_sputum_no_smo), 
  type = 'spearman')

# create node and edge lists
# severe
df_healthy_sputum_no_smo_cor_Pval <- df_healthy_sputum_no_smo_cor$P
df_healthy_sputum_no_smo_COR <- df_healthy_sputum_no_smo_cor$r

# generate edgelist
df_healthy_sputum_no_smo_edges <- melt(
  df_healthy_sputum_no_smo_cor_Pval)
df_healthy_sputum_no_smo_coredges <- melt(
  df_healthy_sputum_no_smo_COR)
df_healthy_sputum_no_smo_edges$COR <- df_healthy_sputum_no_smo_coredges$value
df_healthy_sputum_no_smo_edges$value <- round(
  df_healthy_sputum_no_smo_edges$value, 5)
df_healthy_sputum_no_smo_edges$Label <- row.names(
  df_healthy_sputum_no_smo_coredges)

colnames(df_healthy_sputum_no_smo_edges) <- c(
  'Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_healthy_sputum_no_smo_edges) <- NULL
df_healthy_sputum_no_smo_edges <- df_healthy_sputum_no_smo_edges[
  complete.cases(df_healthy_sputum_no_smo_edges), ]
df_healthy_sputum_no_smo_edges_short <- subset(
  df_healthy_sputum_no_smo_edges, 
  pValue < sig_level & Weight >= weight_val | 
    pValue < sig_level & Weight <= weight_val2)
df_healthy_sputum_no_smo_edges_short$cortype <- ifelse(
  df_healthy_sputum_no_smo_edges_short$Weight > 0, "pos", "neg")
df_healthy_sputum_no_smo_edges_short$Genus <- df_healthy_sputum_no_smo_edges_short$Source
df_healthy_sputum_no_smo_edges_short <- separate(
  df_healthy_sputum_no_smo_edges_short,
  col = 'Genus', into=c('Genus', 'N'), 
  sep='_', remove=TRUE)
df_healthy_sputum_no_smo_edges_short$N <- NULL
rownames(df_healthy_sputum_no_smo_edges_short) <- NULL
df_healthy_sputum_no_smo_edges_short$Weight <- round(
  df_healthy_sputum_no_smo_edges_short$Weight,2)
df_healthy_sputum_no_smo_edges_short$Type <- "Undirected"

df_healthy_sputum_no_smo_nodes <- select(
  df_healthy_sputum_no_smo_edges_short, c(
    "Source", "Genus"))
df_healthy_sputum_no_smo_nodes = df_healthy_sputum_no_smo_nodes[
  !duplicated(df_healthy_sputum_no_smo_nodes$Source),]
df_healthy_sputum_no_smo_nodes <- data.frame(
  df_healthy_sputum_no_smo_nodes)
rownames(df_healthy_sputum_no_smo_nodes) <- NULL
colnames(df_healthy_sputum_no_smo_nodes) <- c(
  "Id", "Label")
df_healthy_sputum_no_smo_nodes1 <- merge(
  df_healthy_sputum_no_smo_nodes, 
  df_healthy_sputum_no_smo_information[, c("Label", "Value")],
  by="Label")

######################################################################

# Output ####
# Files used for network analysis
write.csv(
  df_solid_PI_S_A_edges_short, 
  'output_files/sputum_solid_PI_S_A_edges.csv', 
  row.names = FALSE)
write.csv(
  df_solid_PI_S_A_edges_nodes1, 
  'output_files/sputum_solid_PI_S_A_nodes.csv', 
  row.names = FALSE)
write.csv(
  df_solid_PI_S_B_edges_short, 
  'output_files/sputum_solid_PI_S_B_edges.csv', 
  row.names = FALSE)
write.csv(
  df_solid_PI_S_B_edges_nodes1, 
  'output_files/sputum_solid_PI_S_B_nodes.csv', 
  row.names = FALSE)
write.csv(
  df_solid_PI_S_C_edges_short, 
  'output_files/sputum_solid_PI_S_C_edges.csv', 
  row.names = FALSE)
write.csv(
  df_solid_PI_S_C_edges_nodes1, 
  'output_files/sputum_solid_PI_S_C_nodes.csv', 
  row.names = FALSE)
write.csv(
  df_healthy_pre_edges_short,
  'output_files/coughSwab_illumina_healthy_preschool_edges.csv', 
  row.names = FALSE)
write.csv(
  df_healthy_pre_edges_nodes1, 
  'output_files/coughSwab_illumina_healthy_preschool_nodes.csv', 
  row.names = FALSE)
write.csv(
  df_CF_pre_edges_short,
  'output_files/coughSwab_illumina_CF_preschool_edges.csv', 
  row.names = FALSE)
write.csv(
  df_CF_pre_edges_nodes1, 
  'output_files/coughSwab_illumina_CF_preschool_nodes.csv', 
  row.names = FALSE)
write.csv(
  df_healthy_sputum_no_smo_edges_short,
  'output_files/sputum_solid_healthy_S_C_edges.csv', 
  row.names = FALSE)
write.csv(
  df_healthy_sputum_no_smo_nodes1, 
  'output_files/sputum_solid_healthy_S_C_nodes.csv', 
  row.names = FALSE)

# Files not used for network analysis
#write.csv(
 # df_healthy_counts_cor_edges_short,
  #'output_files/coughSwab_illumina_healthy_C_edges.csv', 
  #row.names = FALSE)
#write.csv(
 # df_healthy_counts_cor_nodes1, 
  #'output_files/coughSwab_illumina_healthy_C_nodes.csv', 
  #row.names = FALSE)
#write.csv(
 # df_water_blank_cor_edges_short, 
  #'output_files/water_blank_edges.csv', 
  #row.names = FALSE)
#write.csv(
 # df_water_blank_cor_nodes1, 
  #'output_files/water_blank_nodes.csv', 
  #row.names = FALSE)
#write.csv(
 # df_swab_blank_cor_edges_short, 
  #'output_files/swab_blank_edges.csv', 
  #row.names = FALSE)
#write.csv(
 # df_swab_blank_cor_nodes1,
  #'output_files/swab_blank_nodes.csv', 
  #row.names = FALSE)
#write.csv(
#  df_SRR_complete_severe_edges_short, 
#  'output_files/sputum_SRR_severe_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_SRR_complete_severe_nodes1, 
 # 'output_files/sputum_SRR_severe_nodes.csv', 
#  row.names = FALSE)
#write.csv(
 # df_SRR_complete_mild_edges_short, 
  #'output_files/sputum_SRR_mild_edges.csv', 
  #row.names = FALSE)
#write.csv(
#  df_SRR_complete_mild_nodes1, 
#  'output_files/sputum_SRR_mild_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_A_edges_short,
#  'output_files/sputum_solid_PS_S_A_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_A_edges_nodes1, 
#  'output_files/sputum_solid_PS_S_A_nodes.csv',
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_B_edges_short, 
#  'output_files/sputum_solid_PS_S_B_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_B_edges_nodes1, 
#  'output_files/sputum_solid_PS_S_B_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_C_edges_short, 
#  'output_files/sputum_solid_PS_S_C_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_S_C_edges_nodes1, 
#  'output_files/sputum_solid_PS_S_C_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_N_A_edges_short,
#  'output_files/nasal_solid_PS_N_A_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PS_N_A_edges_nodes1, 
#  'output_files/nasal_solid_PS_N_A_nodes.csv',
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_A_edges_short, 
#  'output_files/nasal_solid_PI_N_A_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_A_edges_nodes1, 
#  'output_files/nasal_solid_PI_N_A_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_B_edges_short, 
#  'output_files/nasal_solid_PI_N_B_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_B_edges_nodes1, 
#  'output_files/nasal_solid_PI_N_B_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_C_edges_short, 
#  'output_files/nasal_solid_PI_N_C_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_solid_PI_N_C_edges_nodes1, 
#  'output_files/nasal_solid_PI_N_C_nodes.csv', 
#  row.names = FALSE)
#write.csv(
#  df_healthy_sputum_smo_edges_short,
#  'output_files/sputum_solid_smoker_S_C_edges.csv', 
#  row.names = FALSE)
#write.csv(
#  df_healthy_sputum_smo_nodes1, 
#  'output_files/sputum_solid_smoker_S_C_nodes.csv', 
#  row.names = FALSE)


