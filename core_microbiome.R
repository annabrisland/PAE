library(tidyverse)
library(phyloseq)
library(devtools)
#install_github("microbiome/microbiome")
library(microbiome)
library(eulerr)

setwd("~/Desktop/PAE")

## read in files
df = readRDS("PAE_minocycline_filtered_phyloseq.RDS")

ctrl_df <- subset_samples(df, prenatal_group == "Control")
drug_df <- subset_samples(df, prenatal_group == "Ethanol")

#df <- transform_sample_counts(df , function(x) x/sum(x))
table(meta(drug_df)$drug_group, useNA = "always")
df.rel <- microbiome::transform(drug_df, "compositional")

groups <- unique(as.character(meta(df.rel)$drug_group))
print(groups)

list_core <- c() # an empty object to store information

for (n in groups){ # for each variable n in groups
  
  ps.sub <- subset_samples(df.rel, drug_group == n) # Choose sample from group by n
  
  core_m <- core_members(ps.sub,
                         detection = 0.001, # in at least 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each group
  list_core[[n]] <- core_m # add to a list core taxa for each group.
}

mycols <- c(`Adolescent Minocycline`="#d6e2e9", `No drug`="#cbf3f0", `Lactational Minocycline`="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)


core_taxa <- core_members(df, detection=0.001, prevalence = 0.5)

PAE_core_taxa <- prune_taxa(core_taxa, df ) %>% tax_table()
