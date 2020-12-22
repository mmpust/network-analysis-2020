# title: "Network statistics"
# author: "Marie-Madlen Pust"
# date: "22 12 2020"

# clear environment
rm(list=ls())

# load packages
library('readr')
library('rcompanion')
library('ggstatsplot')
library('ggplot2')
library('ggpubr')
library('conover.test')
library('rayshader')

# set global variables
set_base_size = 12

# import gephi output files
healthy_0 <- read_csv(
  "gephi_output/healthy_0_statistics.csv")
healthy_0 <- data.frame(
  healthy_0)
healthy_0$timeset <- NULL
healthy_0$Label <- NULL
healthy_0$Degree_norm <-
  healthy_0$Degree / nrow(healthy_0)
healthy_0$modularity_class <- 
  healthy_0$modularity_class + 1

healthy_C <- read_csv(
  "gephi_output/healthy_C_statistics.csv")
healthy_C <- data.frame(
  healthy_C)
healthy_C$timeset <- NULL
healthy_C$Label <- NULL
healthy_C$Degree_norm <-
  healthy_C$Degree / nrow(healthy_C)
healthy_C$modularity_class <- 
  healthy_C$modularity_class + 1

cf_pi_0 <- read_csv(
  "gephi_output/cf_pi_0_statistics.csv")
cf_pi_0 <- data.frame(
  cf_pi_0)
cf_pi_0$timeset <- NULL
cf_pi_0$Label <- NULL
cf_pi_0$Degree_norm <-
  cf_pi_0$Degree / nrow(cf_pi_0)
cf_pi_0$modularity_class <- 
  cf_pi_0$modularity_class + 1

cf_pi_A <-read_csv(
  "gephi_output/cf_pi_A_statistics.csv")
cf_pi_A <- data.frame(
  cf_pi_A)
cf_pi_A$timeset <- NULL
cf_pi_A$Label <- NULL
cf_pi_A$Degree_norm <-
  cf_pi_A$Degree / nrow(cf_pi_A)
cf_pi_A$modularity_class <- 
  cf_pi_A$modularity_class + 1

cf_pi_B <-read_csv(
  "gephi_output/cf_pi_B_statistics.csv")
cf_pi_B <- data.frame(
  cf_pi_B)
cf_pi_B$timeset <- NULL
cf_pi_B$Label <- NULL
cf_pi_B$Degree_norm <-
  cf_pi_B$Degree / nrow(cf_pi_B)
cf_pi_B$modularity_class <- 
  cf_pi_B$modularity_class + 1

cf_pi_C <- read_csv(
  "gephi_output/cf_pi_C_statistics.csv")
cf_pi_C <- data.frame(
  cf_pi_C)
cf_pi_C$timeset <- NULL
cf_pi_C$Label <- NULL
cf_pi_C$Degree_norm <-
  cf_pi_C$Degree / nrow(cf_pi_C)
cf_pi_C$modularity_class <- 
  cf_pi_C$modularity_class + 1

# merge all tables
healthy_0$ageGroup <- "4-6 (H)"
healthy_0$state <- "healthy"
healthy_C$ageGroup <- ">17 (H)"
healthy_C$state <- "healthy"
healthy_C$Weighted.Degree <- NULL
cf_pi_0$ageGroup <- "4-6 (CF)"
cf_pi_0$state <- "cf"
cf_pi_0$Weighted.Degree <- NULL
cf_pi_A$ageGroup <- "7-12 (CF)"
cf_pi_A$state <- "cf"
cf_pi_A$Weighted.Degree <- NULL
cf_pi_B$ageGroup <- "13-17 (CF)"
cf_pi_B$state <- "cf"
cf_pi_B$Weighted.Degree <- NULL
cf_pi_C$ageGroup <- ">17 (CF)"
cf_pi_C$state <- "cf"
cf_pi_C$Weighted.Degree <- NULL

table_all <- (data.frame(
  rbind(healthy_0, healthy_C,
        cf_pi_0, cf_pi_A, 
        cf_pi_B, cf_pi_C)))
table_all$ageGroup <- factor(as.character(
  table_all$ageGroup), 
  levels = c(
    "4-6 (H)", "4-6 (CF)", 
    "7-12 (CF)", "13-17 (CF)", 
    ">17 (CF)", ">17 (H)"))


# statistics
# degree centrality ####
degree_plot <- ggplot(table_all, aes(
  x=ageGroup,
  y=Degree_norm,
  color=state)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.05,
    alpha = 0.4) +
  theme_pubr(
    border=TRUE,
    base_size = set_base_size,
    legend = "none") +
  xlab(
    " ") +
  ylab(
    "Degree\n") +
  ylim(0.0,1.0) +
  scale_x_discrete(
    labels=c(
      "4-6", "4-6", 
      "7-15", "16-23", 
      ">23", ">23")) +
  theme(
    axis.text.x = element_text(angle = 90)) +
  scale_color_manual(
    values=c("darkred", 
             "darkblue"))

kruskal.test(
  table_all$Degree, 
  table_all$ageGroup)

epsilonSquared(
  table_all$Degree,
  table_all$ageGroup, 
  ci = TRUE)

conover.test(
  table_all$Degree,
  g=table_all$ageGroup,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

# statistics
# modularity class ####
modularity_plot <- ggplot(table_all, aes(
  x=ageGroup,
  y=modularity_class,
  color=state)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.05,
    alpha = 0.4) +
  theme_pubr(
    border=TRUE,
    base_size = set_base_size,
    legend = "none") +
  xlab(
    " ") +
  ylab(
    "Modularity class\n") +
  scale_x_discrete(
    labels=c(
      "4-6", "4-6", 
      "7-15", "16-23", 
      ">23", ">23")) +
  theme(
    axis.text.x = element_text(angle = 90)) +
  scale_color_manual(
    values=c("darkred", 
             "darkblue"))

kruskal.test(
  table_all$modularity_class, 
  table_all$ageGroup)

epsilonSquared(
  table_all$modularity_class,
  table_all$ageGroup, 
  ci = TRUE)

conover.test(
  table_all$modularity_class,
  g=table_all$ageGroup,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

# statistics
# eigencentrality ####
eigenvector_plot <- ggplot(table_all, aes(
  x=ageGroup,
  y=eigencentrality,
  color=state)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.05,
    alpha = 0.4) +
  theme_pubr(
    border=TRUE,
    base_size = set_base_size,
    legend = "none") +
  xlab(
    "Age (in years)") +
  ylab(
    "Eigencentrality\n") +
  ylim(0.0,1.0) +
  scale_x_discrete(
    labels=c(
      "4-6", "4-6", 
      "7-15", "16-23", 
      ">23", ">23")) +
  theme(
    axis.text.x = element_text(angle = 90)) +
  scale_color_manual(
    values=c("darkred", 
             "darkblue"))

kruskal.test(
  table_all$eigencentrality, 
  table_all$ageGroup)

epsilonSquared(
  table_all$eigencentrality,
  table_all$ageGroup, 
  ci = TRUE)

conover.test(
  table_all$eigencentrality,
  g=table_all$ageGroup,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.05, 
  altp=TRUE)


# differences in eigencentrality 
# across modularity classes
# cf-0 ####
cf_pi_0$modularity_class <- as.factor(as.numeric(
  cf_pi_0$modularity_class))
levels(cf_pi_0$modularity_class)
kruskal.test(
  cf_pi_0$eigencentrality, 
  cf_pi_0$modularity_class)

epsilonSquared(
  cf_pi_0$eigencentrality,
  cf_pi_0$modularity_class, 
  ci = TRUE)

conover.test(
  cf_pi_0$eigencentrality,
  g=cf_pi_0$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.05, 
  altp=TRUE)

comparision_cf_0_mod <- list(
  c("4", "3"),
  c("4", "2"),
  c("4", "1"),
  c("3", "2"),
  c("3", "1"),
  c("2", "1"))

cf_0_mod_eig_plot <- ggplot(
  data=cf_pi_0, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_boxplot() +
  geom_jitter(aes(
    colour=value,
    size=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    "Modularity class") +
  ylab(
    "Eigencentrality") +
  ylim(
    0.0, 2.3) +
  scale_colour_gradient(
    low="orange",
    high="darkred",
    labels=scales::number_format(accuracy = 1)) +
  scale_size(
    labels=scales::number_format(accuracy = 1)) +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance")) +
  stat_compare_means(
    comparisons=comparision_cf_0_mod,
    label.y = c(
      1.2, 1.4, 
      1.6, 1.8, 
      2.0, 2.2),
    label = "p.signif") 

# correlation between abundance and modularity class
kruskal.test(
  cf_pi_0$value,
  as.factor(as.numeric(
    cf_pi_0$modularity_class)))
epsilonSquared(
  cf_pi_0$value,
  as.factor(as.numeric(
    cf_pi_0$modularity_class)),
  ci = TRUE)

# crrelation between abundance and eigencentrality
cor.test(
  cf_pi_0$value, 
  cf_pi_0$eigencentrality, 
  method = "spearman")
# no correlation between absolute abundance and eigencentrality.

# correlation between abundance and degree centrality
cor.test(
  cf_pi_0$value, 
  cf_pi_0$Degree_norm, 
  method = "spearman")

# cf-A ####
cf_pi_A$modularity_class <- as.factor(as.numeric(
  cf_pi_A$modularity_class))
levels(cf_pi_A$modularity_class)
kruskal.test(
  cf_pi_A$eigencentrality, 
  cf_pi_A$modularity_class)

epsilonSquared(
  cf_pi_A$eigencentrality,
  cf_pi_A$modularity_class, 
  ci = TRUE)

conover.test(
  cf_pi_A$eigencentrality,
  g=cf_pi_A$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

comparision_cf_A_mod <- list(
  c("1", "2"),
  c("1", "4"),
  c("1", "5"),
  c("3", "4"),
  c("3", "5"),
  c("4", "5"))

cf_pi_A$value <-
  round(cf_pi_A$value,1)
cf_A_mod_eig_plot <- ggplot(
  data=cf_pi_A, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_boxplot(
  ) +
  geom_jitter(aes(
    colour=value,
    size=value,
    fill=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    "Modularity class") +
  ylab(
    "Eigencentrality") +
  ylim(
    0.0, 2.3) +
  scale_colour_gradient(
    low="orange",
    high="darkred") +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance"),
    fill = guide_legend("Abundance")) +
  stat_compare_means(
    comparisons=comparision_cf_A_mod,
    label.y = c(
      1.2, 1.4, 1.6, 1.8, 
      2.0, 2.2),
    label = "p.signif") 

#  correlation between abundance and modularity class
kruskal.test(
  cf_pi_A$value,
  as.factor(as.numeric(
    cf_pi_A$modularity_class)))


# cf-B ####
cf_pi_B$modularity_class <- as.factor(as.numeric(
  cf_pi_B$modularity_class))
levels(cf_pi_B$modularity_class)
kruskal.test(
  cf_pi_B$eigencentrality, 
  cf_pi_B$modularity_class)

#epsilonSquared(
 # cf_pi_B$eigencentrality,
  #cf_pi_B$modularity_class, 
  #ci = TRUE)

#conover.test(
#  cf_pi_B$eigencentrality,
#  g=cf_pi_B$modularity_class,
#  method="bh", 
#  kw=TRUE, 
#  label=TRUE, 
#  wrap=TRUE, 
#  table=TRUE, 
#  list=FALSE, 
#  rmc=FALSE, 
#  alpha=0.01, 
#  altp=TRUE)

cf_pi_B$value <-
  round(cf_pi_B$value,1)
cf_B_mod_eig_plot <- ggplot(
  data=cf_pi_B, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_boxplot(
  ) +
  geom_jitter(aes(
    colour=value,
    size=value,
    fill=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    " ") +
  ylab(
    "Eigencentrality") +
  ylim(
    0.0, 2.3) +
  scale_colour_gradient(
    low="orange",
    high="darkred") +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance"),
    fill = guide_legend("Abundance")) 

# correlation between abundance and eigencentrality
#cor.test(
#  cf_pi_B$value, 
#  cf_pi_B$eigencentrality, 
#  method = "spearman")
# no correlation between absolute abundance and eigencentrality.
# p-value ~ 0.15, rho = -0.46, S = 321.8


# cf-C ####
cf_pi_C$modularity_class <- as.factor(as.numeric(
  cf_pi_C$modularity_class))
levels(cf_pi_C$modularity_class)
kruskal.test(
  cf_pi_C$eigencentrality, 
  cf_pi_C$modularity_class)

epsilonSquared(
  cf_pi_C$eigencentrality,
  cf_pi_C$modularity_class, 
  ci = TRUE)

conover.test(
  cf_pi_C$eigencentrality,
  g=cf_pi_C$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.05, 
  altp=TRUE)

comparision_cf_C_mod <- list(
  c("1", "2"),
  c("1", "3"),
  c("1", "4"),
  c("2", "3"))

cf_pi_C$value <-
  round(cf_pi_C$value,1)
cf_C_mod_eig_plot <- ggplot(
  data=cf_pi_C, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_boxplot(
  ) +
  geom_jitter(aes(
    colour=value,
    size=value,
    fill=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    "Modularity class") +
  ylab(
    " ") +
  ylim(
    0.0, 2.3) +
  scale_colour_gradient(
    low="orange",
    high="darkred") +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance"),
    fill = guide_legend("Abundance")) +
  stat_compare_means(
    comparisons=comparision_cf_C_mod,
    label.y = c(
      1.4, 1.6, 
      1.8, 2.0),
    label = "p.signif") 

# correlation between abundance and modularity class
kruskal.test(
  cf_pi_C$value,
  as.factor(as.numeric(
    cf_pi_C$modularity_class)))

epsilonSquared(
  cf_pi_C$value,
  as.factor(as.numeric(
    cf_pi_C$modularity_class)),
  ci = TRUE)

# correlation between abundance and eigencentrality
cor.test(
  cf_pi_C$value, 
  cf_pi_C$eigencentrality, 
  method = "spearman")

# correlation between abundance and degree centrality
cor.test(
  cf_pi_C$value, 
  cf_pi_C$Degree_norm, 
  method = "spearman")


# healthy_0 ####
healthy_0$modularity_class <- as.factor(as.numeric(
  healthy_0$modularity_class))
levels(healthy_0$modularity_class)
kruskal.test(
  healthy_0$eigencentrality, 
  healthy_0$modularity_class)

epsilonSquared(
  healthy_0$eigencentrality,
  healthy_0$modularity_class, 
  ci = TRUE)

conover.test(
  healthy_0$eigencentrality,
  g=healthy_0$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  list=FALSE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

#comparision_healthy_0_mod <- list(
#  c("1", "2"),
#  c("1", "3"),
#  c("1", "4"),
#  c("2", "4"),
#  c("3", "4"))

healthy_0$value <-
  round(healthy_0$value,1)

H_0_mod_eig_plot <- ggplot(
  data=healthy_0, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_boxplot(
  ) +
  geom_jitter(aes(
    colour=value,
    size=value,
    fill=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    " ") +
  ylab(
    "Eigencentrality") +
  ylim(
    0.0, 2.2) +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance"),
    fill = guide_legend("Abundance")) #+
 # stat_compare_means(
  #  comparisons=comparision_healthy_0_mod,
   # label.y = c(
    #  1.2, 1.4, 1.6, 
     # 1.8, 2.0),
   # label = "p.signif") 

# correlation between abundance and eigencentrality
cor.test(
  healthy_0$value, 
  healthy_0$eigencentrality, 
  method = "spearman")
# no correlation between absolute abundance and eigencentrality.
# p-value ~ 0.8, rho = -0.007, S = 114125428

# correlation between abundance and degree centrality
cor.test(
  healthy_0$value, 
  healthy_0$Degree_norm, 
  method = "spearman")
# no correlation between absolute abundance and eigencentrality.
# p-value ~ 0.8, rho = -0.007, S = 114125428
# 
# correlation between abundance and modularity class
kruskal.test(
  healthy_0$value, 
  healthy_0$modularity_class)

epsilonSquared(
  healthy_0$value,
  healthy_0$modularity_class, 
  ci = TRUE)

# healthy_C ####
healthy_C$modularity_class <- as.factor(as.numeric(
  healthy_C$modularity_class))
levels(healthy_C$modularity_class)
kruskal.test(
  healthy_C$eigencentrality, 
  healthy_C$modularity_class)

epsilonSquared(
  healthy_C$eigencentrality,
  healthy_C$modularity_class, 
  ci = TRUE)

conover.test(
  healthy_C$eigencentrality,
  g=healthy_C$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=TRUE, 
  table=TRUE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

comparision_healthy_C_mod <- list(
  c("1", "2"),
  c("1", "3"),
  c("1", "4"),
  c("1", "5"),
  c("1", "6"),
  c("1", "7"),
  c("1", "8"))

healthy_C$value <-
  round(healthy_C$value,1)
H_C_mod_eig_plot <- ggplot(
  data=healthy_C, aes(
    x=modularity_class,
    y=eigencentrality)) +
  geom_smooth(
  ) +
  geom_jitter(aes(
    colour=value,
    size=value,
    fill=value),
    width = 0.1, 
    alpha=0.5) +
  theme_pubr(
    border=TRUE,
    legend="right",
    base_size = set_base_size) +
  theme(
    axis.text.x = element_blank()) +
  xlab(
    " ") +
  ylab(
    " ") +
  ylim(
    0.0, 2.2) +
  guides(
    colour=guide_legend("Abundance"), 
    size = guide_legend("Abundance"),
    fill = guide_legend("Abundance")) 

# correlation between abundance and eigencentrality
cor.test(
  healthy_C$value, 
  healthy_C$eigencentrality, 
  method = "spearman")
# no correlation between absolute abundance and eigencentrality.
# p-value ~ 0.3, rho = -0.03, S = 81504520

# correlation between abundance and degree centrality
cor.test(
  healthy_C$value, 
  healthy_C$Degree_norm, 
  method = "spearman")
# no correlation between absolute abundance and eigencentrality.
# p-value ~ 0.2, rho = -0.04, S = 81906173
# 
# correlation between abundance and modularity class
kruskal.test(
  healthy_C$value, 
  healthy_C$modularity_class)

epsilonSquared(
  healthy_C$value,
  healthy_C$modularity_class, 
  ci = TRUE)

test <-
  conover.test(
  healthy_C$value,
  g=healthy_C$modularity_class,
  method="bh", 
  kw=TRUE, 
  label=TRUE, 
  wrap=FALSE, 
  table=TRUE, 
  rmc=FALSE, 
  alpha=0.01, 
  altp=TRUE)

# save images to working directory
network_stats1 <- ggarrange(
  degree_plot, 
  eigenvector_plot, 
  modularity_plot,
  nrow = 1, 
  labels = c("A", "B", "C"))
network_stats1


network_stats2 <- ggarrange(
  H_0_mod_eig_plot,
  H_C_mod_eig_plot,
  cf_0_mod_eig_plot,
  cf_C_mod_eig_plot, nrow=2, ncol=2,
  labels = c("A", "B", "C", "D"))
network_stats2

ggsave("Figure_02.tif",
       network_stats1,
       device="tiff",
       scale=1,
       width = 7.57,
       height = 4.07,
       dpi=600)

ggsave("Figure_03.tif",
       network_stats2,
       device="tiff",
       scale=1,
       width = 7.57,
       height = 5.34,
       dpi=600)
