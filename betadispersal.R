library(phyloseq)
library(tidyverse)
library(vegan)
library(devtools)

set.seed(333)

## read in files
PAE = readRDS("PAE_minocycline_rarefied_df.RDS")

project_bray.rarefied <- phyloseq::distance(PAE, method = "bray")
sample_df <- data.frame(sample_data(PAE))
beta.FACTOR2 <- betadisper(project_bray.rarefied, sample_df$drug_group) 
b1.2 = permutest(beta.FACTOR2) 
b1.2
