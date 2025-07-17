# set up environment
  setwd("/Users/leenakwak/Desktop/R projects/SDoH Paper")

# install packages
  install.packages("devtools")
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


# load packages -----
  library(vegan)
  library(pairwiseAdonis)
  library(readr)
  library(readxl)
  library(tidyverse)
  library(dplyr)    
  library(ggplot2)
  library(ggtext)
  library(ggpubr)
  library(pairwiseAdonis)

# read in data (converted .qza files to .tsv manually) -----

  weighted_unifrac = as.dist(read.table("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/unweighted-unifrac-distance-matrix.tsv", header = T))
  unweighted_unifrac = as.dist(read.table("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/weighted-unifrac-distance-matrix.tsv", header = T))
  bray_3300 = as.dist(read.table("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/bray-curtis-distance-matrix.tsv", header = T))
  jaccard_3300 = as.dist(read.table("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/jaccard-distance-matrix.tsv", header = T))
  
  shannon <- read_tsv("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/shannon-diversity.tsv")
  pielou_evenness <- read_tsv("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/pielou-evenness.tsv")
  colnames(pielou_evenness)[1] <- "sample_id"
  faith_pd <- read_tsv("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/faith-pd.tsv")
  observed_features <- read_tsv("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/observed-features.tsv")
  colnames(observed_features)[1] <- "sample_id"
  
  glm_phyla <- read.csv("/Users/leenakwak/Desktop/R projects/SDoH Paper/data/glm_phyla.csv")
  
  metadata <- read_excel(path="/Users/leenakwak/Desktop/R projects/SDoH Paper/data/REACH_metadata_5.xlsx")
  
  set.seed(1018)
  
  
# income + beta diversity -----
  
  mds_unweighted_unifrac<-metaMDS(unweighted_unifrac, k=2, trymax=1000)
  mds_unweighted_unifrac_points<-mds_unweighted_unifrac$points
  mds_unweighted_unifrac_points2<-merge(x=mds_unweighted_unifrac_points, y = metadata, 
                               by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical income_range to descriptive labels
  income_labels <- c("0" = "Less than $10,000", 
                     "1" = "$10,000 - $19,999", 
                     "2" = "$20,000 - $34,999", 
                     "3" = "$35,000 - $49,999", 
                     "4" = "$50,000 - $74,999", 
                     "5" = "$75,000 - $99,999", 
                     "6" = "$100,000+")
  
  mds_unweighted_unifrac_points2$Income_range <- factor(mds_unweighted_unifrac_points2$Income_range, 
                                                        levels = names(income_labels), 
                                                        labels = income_labels)
  
  unweighted_unifrac_plot_income_range <- ggplot(data = subset(mds_unweighted_unifrac_points2, !is.na(Income_range)), 
                                                 aes(x = MDS1, y = MDS2, color = Income_range, shape = Income_range)) +
    labs(shape = "Income Range", color = "Income Range") + 
    scale_color_manual(values = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#000000")) +
    scale_shape_manual(values = c(0, 2, 3, 8, 5, 6)) +
    geom_point(size = 3) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(unweighted_unifrac_plot_income_range)
  ggsave("income range unweighted unifrac.jpeg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# income + alpha diversity-----
  colnames(shannon)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = shannon, by.x = "Sample_ID", by.y = "sample_id")

  
# convert numerical Income_range to descriptive labels
  income_labels <- c("0" = "Less than $10,000", 
                     "1" = "$10,000 - $19,999", 
                     "2" = "$20,000 - $34,999", 
                     "3" = "$35,000 - $49,999", 
                     "4" = "$50,000 - $74,999", 
                     "5" = "$75,000 - $99,999", 
                     "6" = "$100,000+")
  
  
  alpha$Income_range <- factor(alpha$Income_range, 
                               levels = names(income_labels), 
                               labels = income_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(Income_range)),  
                     x = "Income_range", 
                     y = "shannon_entropy.y", 
                     color = "Income_range", 
                     palette = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#000000"), 
                     add = "jitter", 
                     ylab = "Shannon Diversity") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Income Range")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))

  print(plot3)
  ggsave("income range shannon.jpg", width = 6, height = 4, device='jpeg', dpi=300)
  
# education + beta diversity-----
  
  mds_jaccard_3300<-metaMDS(jaccard_3300, k=2, trymax=1000)
  mds_jaccard_3300_points<-mds_jaccard_3300$points
  mds_jaccard_3300_points2<-merge(x=mds_jaccard_3300_points, y = metadata, 
                                        by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical Education to descriptive labels
  education_labels <- c("1" = "Less than high school", 
                        "2" = "High school", 
                        "3" = "College, no degree", 
                        "4" = "Associate degree", 
                        "5" = "Bachelor's degree", 
                        "6" = "Master's degree")
  
  mds_jaccard_3300_points2$Education <- factor(mds_jaccard_3300_points2$Education, 
                                                        levels = names(education_labels), 
                                                        labels = education_labels)
  
  jaccard_3300_plot_education <- ggplot(data = subset(mds_jaccard_3300_points2, !is.na(Education)), 
                                                 aes(x = MDS1, y = MDS2, color = Education, shape = Education)) +
    labs(shape = "Education", color = "Education") + 
    scale_color_manual(values = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#000000")) +
    scale_shape_manual(values = c(0, 2, 3, 8, 5, 6)) +
    geom_point(size = 3) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(jaccard_3300_plot_education)
  ggsave("education jaccard.jpeg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# education + alpha diversity -----
  colnames(observed_features)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = observed_features, by.x = "Sample_ID", by.y = "sample_id")
  
  
# convert numerical Education to descriptive labels
  education_labels <- c("1" = "Less than high school", 
                     "2" = "High school", 
                     "3" = "College, no degree", 
                     "4" = "Associate degree", 
                     "5" = "Bachelor's degree", 
                     "6" = "Master's degree")
  
  
  alpha$Education <- factor(alpha$Education, 
                               levels = names(education_labels), 
                               labels = education_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(Education)),  
                     x = "Education", 
                     y = "observed_features.y", 
                     color = "Education", 
                     palette = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#000000"), 
                     add = "jitter", 
                     ylab = "Observed Features") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Education")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(plot3)
  ggsave("education observed features.jpg", width = 6, height = 4, device='jpeg', dpi=300)

  
# healthcare access + beta diversity-----
  
  mds_weighted_unifrac<-metaMDS(weighted_unifrac, k=2, trymax=1000)
  mds_weighted_unifrac_points<-mds_weighted_unifrac$points
  mds_weighted_unifrac_points2<-merge(x=mds_weighted_unifrac_points, y = metadata, 
                                  by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical healthcare_access to descriptive labels
  healthcare_labels <- c("0" = "Sufficient", 
                         "1" = "Limited")
  
  mds_weighted_unifrac_points2$healthcare_access <- factor(mds_weighted_unifrac_points2$healthcare_access, 
                                               levels = names(healthcare_labels), 
                                               labels = healthcare_labels)
  
  weighted_unifrac_plot_healthcare <- ggplot(data = subset(mds_weighted_unifrac_points2, !is.na(healthcare_access)), 
                                        aes(x = MDS1, y = MDS2, color = healthcare_access, shape = healthcare_access)) +
    labs(shape = "Healthcare Access", color = "Healthcare Access") + 
    scale_color_manual(values = c("#0072B2","#CC79A7","#E69F00","#D55E00","#56B4E9", "#009E73", "#000000")) +
    scale_shape_manual(values = c(16, 17)) +
    geom_point(size = 3) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1))) +
    stat_ellipse(aes(group = healthcare_access), type = "t", level = 0.9)
  
  print(weighted_unifrac_plot_healthcare)
  ggsave("healthcare weighted unifrac.jpeg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# healthcare access + alpha diversity-----
  colnames(observed_features)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = observed_features, by.x = "Sample_ID", by.y = "sample_id")
  
# convert numerical healthcare_access to descriptive labels
  healthcare_labels <- c("0" = "Sufficient", 
                        "1" = "Limited")
  
  alpha$healthcare_access <- factor(alpha$healthcare_access, 
                            levels = names(healthcare_labels), 
                            labels = healthcare_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(healthcare_access)),  
                     x = "healthcare_access", 
                     y = "observed_features.y", 
                     color = "healthcare_access", 
                     palette = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#000000"), 
                     add = "jitter", 
                     ylab = "Observed Features") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Healthcare Access")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(plot3)
  ggsave("observed features healthcare access.jpg", width = 6, height = 4, device='jpeg', dpi=300)
  

# food insecurity + beta diversity -----
  
  mds_bray_3300<-metaMDS(bray_3300, k=2, trymax=1000)
  mds_bray_3300_points<-mds_bray_3300$points
  mds_bray_3300_points2<-merge(x=mds_bray_3300_points, y = metadata, 
                                        by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical food_insecurity to descriptive labels
  food_insecurity_labels <- c("1" = "Did not skip/reduce meals", 
                              "2" = "Skipped/reduced 1-2 months",
                              "3" = "Skipped/reduced 3-6 months",
                              "4" = "Skipped/reduced 7-9 months",
                              "5" = "Skipped/reduced 10-12 months")
  
  mds_bray_3300_points2$food_insecurity <- factor(mds_bray_3300_points2$food_insecurity, 
                                                        levels = names(food_insecurity_labels), 
                                                        labels = food_insecurity_labels)
  
  bray_3300_plot_food_insecurity <- ggplot(data = subset(mds_bray_3300_points2, !is.na(food_insecurity)), 
                                                 aes(x = MDS1, y = MDS2, color = food_insecurity, shape = food_insecurity)) +
    labs(shape = "Food Insecurity", color = "Food Insecurity") + 
    scale_color_manual(values = c("#E69F00","#CC79A7","#009E73", "#0072B2","#D55E00","#56B4E9", "#000000")) +
    scale_shape_manual(values = c(0, 2, 3, 8, 5, 6)) +
    geom_point(size = 3) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1))) +
    stat_ellipse(aes(group = food_insecurity), type = "t", level = 0.9)
  
  print(bray_3300_plot_food_insecurity)
  ggsave("food insecurity bray curtis.jpeg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# food insecurity + alpha diversity-----
  colnames(observed_features)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = observed_features, by.x = "Sample_ID", by.y = "sample_id")
  
  
  # convert numerical food_insecurity to descriptive labels
  food_insecurity_labels <- c("1" = "Does not skip meals", 
                         "2" = "1-2 months",
                         "3" = "3-6 months",
                         "4" = "7-9 months",
                         "5" = "10-12 months")
  
  
  alpha$food_insecurity <- factor(alpha$food_insecurity, 
                                    levels = names(food_insecurity_labels), 
                                    labels = food_insecurity_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(food_insecurity)),  
                     x = "food_insecurity", 
                     y = "observed_features.y", 
                     color = "food_insecurity", 
                     palette = c("#E69F00","#CC79A7","#56B4E9", "#009E73", "#0072B2", "#000000"), 
                     add = "jitter", 
                     ylab = "Observed Features") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Food Security")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(plot3)
  ggsave("observed features food security.jpg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# overcrowding + beta diversity -----
  
  mds_bray_3300<-metaMDS(bray_3300, k=2, trymax=1000)
  mds_bray_3300_points<-mds_bray_3300$points
  mds_bray_3300_points2<-merge(x=mds_bray_3300_points, y = metadata, 
                                      by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical overcrowding_score to descriptive labels
  overcrowding_labels <- c("1" = "No", 
                           "2" = "Yes")
  
  mds_bray_3300_points2$overcrowding_score <- factor(mds_bray_3300_points2$overcrowding_score, 
                                                           levels = names(overcrowding_labels), 
                                                           labels = overcrowding_labels)
  
  bray_3300_plot_overcrowded <- ggplot(data = subset(mds_bray_3300_points2, !is.na(overcrowding_score)), 
                                             aes(x = MDS1, y = MDS2, color = overcrowding_score, shape = overcrowding_score)) +
    labs(shape = "Overcrowded", color = "Overcrowded") + 
    scale_color_manual(values = c("#0072B2","#CC79A7","#E69F00","#D55E00","#56B4E9", "#009E73", "#000000")) +
    scale_shape_manual(values = c(16, 17)) +
    geom_point(size = 3) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1))) +
    stat_ellipse(aes(group = overcrowding_score), type = "t", level = 0.9)
  
  print(bray_3300_plot_overcrowded)
  ggsave("overcrowded bray curtis.jpeg", width = 6, height = 4, device='jpeg', dpi=300)
  
    
# overcrowding + alpha diversity -----
  colnames(observed_features)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = observed_features, by.x = "Sample_ID", by.y = "sample_id")
  
  
# convert numerical overcrowding_score to descriptive labels
  overcrowding_labels <- c("1" = "No", 
                         "2" = "Yes")
  
  
  alpha$overcrowding_score <- factor(alpha$overcrowding_score, 
                                    levels = names(overcrowding_labels), 
                                    labels = overcrowding_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(overcrowding_score)),  
                     x = "overcrowding_score", 
                     y = "observed_features.y", 
                     color = "overcrowding_score", 
                     palette = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#000000"), 
                     add = "jitter", 
                     ylab = "Observed Features") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Overcrowded")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(plot3)
  ggsave("observed features overcrowding.jpg", width = 6, height = 4, device='jpeg', dpi=300)
  
  
# employment + beta diversity -----
  
  mds_bray_3300<-metaMDS(bray_3300, k=2, trymax=1000)
  mds_bray_3300_points<-mds_bray_3300$points
  mds_bray_3300_points2<-merge(x=mds_bray_3300_points, y = metadata, 
                                        by.x = "row.names", by.y = "Sample_ID")
  
# convert numerical Employment_status to descriptive labels
  employment_labels <- c("0" = "Employed full time", 
                     "1" = "Employed part time", 
                     "2" = "Unemployed, looking for work", 
                     "3" = "Unemployed, not looking for work", 
                     "4" = "Student", 
                     "5" = "Retired", 
                     "6" = "Stay-at-home parent",
                     "7" = "Self-emplyed",
                     "8" = "Unable to work")
  
  mds_bray_3300_points2$Employment_status <- factor(mds_bray_3300_points2$Employment_status, 
                                                        levels = names(employment_labels), 
                                                        labels = employment_labels)
  
  bray_3300_plot_employment <- ggplot(data = subset(mds_bray_3300_points2, !is.na(Employment_status)), 
                                                 aes(x = MDS1, y = MDS2, color = Employment_status, shape = Employment_status)) +
    labs(shape = "Employment Status", color = "Employment Status") + 
    scale_color_manual(values = c("#E69F00","#CC79A7","#D55E00","#56B4E9", "#009E73", "#0072B2", "#CCB974", "#F0E442", "#999999")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 8, 5, 6, 7)) +
    geom_point(size = 2) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          legend.key = element_blank()) +
    theme(axis.title.x = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2)),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(bray_3300_plot_employment)
  ggsave("employment bray.jpeg", width = 6.5, height = 4, device='jpeg', dpi=300)
  
  
# employment + alpha diversity-----
  colnames(observed_features)[1] <- "sample_id"
  alpha <- merge(x = metadata, y = observed_features, by.x = "Sample_ID", by.y = "sample_id")
  
  
# convert numerical Employment_status to descriptive labels
  employment_labels <- c("0" = "Employed full time", 
                         "1" = "Employed part time", 
                         "2" = "Unemployed, looking for work", 
                         "3" = "Unemployed, not looking for work", 
                         "4" = "Student", 
                         "5" = "Retired", 
                         "6" = "Stay-at-home parent",
                         "7" = "Self-emplyed",
                         "8" = "Unable to work")
  
  
  alpha$Employment_status <- factor(alpha$Employment_status, 
                               levels = names(employment_labels), 
                               labels = employment_labels)
  
  plot3 <- ggboxplot(data = subset(alpha, !is.na(Employment_status)),  
                     x = "Employment_status", 
                     y = "observed_features.y", 
                     color = "Employment_status", 
                     palette = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00","#F0E442", "#999999", "#000000"), 
                     add = "jitter", 
                     ylab = "Observed Features") + 
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5)) +  
    guides(color = guide_legend(title = "Employment Status")) + 
    theme(text = element_text(family = "Tahoma")) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          plot.title = element_text(size = rel(2), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)))
  
  print(plot3)
  ggsave("employment observed features.jpg", width = 6.5, height = 4, device='jpeg', dpi=300)
  
  
  