library(phyloseq)
library(tidyverse)
library(vegan)
#library(ggrare)
library(stats)

## tell R where to get the files
setwd("~/Desktop/PAE")

## read in files
df = readRDS("pae_minocycline_unfiltered_phyloseq.RDS")
df

##### FILTERING - REMOVE OFF-TAGET TAXA AND MAJOR CONTAMINATION FROM TAXONOMY FILE #####
df = subset_taxa(df, Kingdom="Bacteria")

##### FILTERING - REMOVE SAMPLES WITH LESS THAN 1000 READS #####
sample_sums(df)
plot(sort(sample_sums(df)))

## add number of total sequences in a sample (Read depth) to the metadata 
df@sam_data$read_depth_noofftargets = sample_sums(df) 

## check which samples have less than 1000 reads
which(df@sam_data$read_depth_noofftargets < 1000)

## remove the samples by your threshold 
df.pruned <- prune_samples(sample_sums(df) >= 1000, df)

## write file to know which samples were lost here. This is important for the methods section. 
df.below1000 <- prune_samples(sample_sums(df) < 1000, df)
df.below1000 = as.matrix(df.below1000@sam_data)
write_rds(df.below1000, "SF_samples_less_than_1000.csv")

##### FILTERING - REMOVE INDIVIDUAL ASVS WITH LESS THAN 100 READS #####
## extract OTU dataframe from phyloseq object
otu.pruned <- as.data.frame(as.matrix(otu_table(df.pruned)))

## remove ASVs (rows) with less than 100 reads across whole dataset but keep all samples
otu.pruned$rowsum = rowSums(otu.pruned)

## remove low frequency asvs
otu.pruned = subset(otu.pruned, otu.pruned$rowsum>100)

## remove rowsum column from your OTU table 
otu = subset(otu.pruned, select=-c(rowsum))


##### FILTERING - REMOVE ASVs FOUND IN 2 SAMPLES OR LESS #####

# has sample ID as column name and asv is as row name. Needs to be this way to use richness function 
# function to calculate richness, sums along a row (OTU)
richness = function(x){return(sum(x>0))}

## calculate richness on entire dataframe
otu$richness = apply(otu,1,richness) # use all columns of otu dataframe
summary(otu$richness)

## remove OTU (rows) with richness 2 or less
otu = subset(otu, otu$richness> 2)
## check that it worked
summary(otu$richness)
## remove richness column
otu = subset(otu, select=-c(richness))


##### DENOISING - MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES 5 OR LESS 0 #####
otu <- mutate_all(otu, funs(ifelse(. < 5, 0, .)))

##### CREATE AND READ BACK IN FILTERED BUT NOT RAREFIED PHYLOSEQ OBJECT ######

## format for phyloseq
otu_mat = as.matrix(otu)
## you need to check if the rows are the ASV sequence (TRUE) or not (FALSE)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

## use the taxonomy table you had in your original phyloseq object
tax_mat = df.pruned@tax_table
TAX = tax_table(tax_mat)  

## get your metadata 
samples = as.data.frame(as.matrix(df.pruned@sam_data))

## get metadata ready for phyloseq
samples = sample_data(samples)

# make new phyloseq object
df.prerarefaction = phyloseq(OTU, TAX, samples)
# check that the phyloseq object was made correctly
df.prerarefaction
rank_names(df.prerarefaction)
sample_variables(df.prerarefaction)


##### FINAL DENOISING AND SAVE FILTERED DATA #####

## remove taxa from taxonomy table with empty otus (from denoising)
df.prerarefaction = prune_taxa(taxa_sums(df.prerarefaction)>0, df.prerarefaction)

## get final read depth
df.prerarefaction@sam_data$read_depth_filtered = sample_sums(df.prerarefaction)

histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)

## save your filtered but not rarefied phyloseq object 
write_rds(df.prerarefaction, "PAE_minocycline_filtered_phyloseq.RDS")


##### SIMPLE RAREFACTION - NORMALIZING THE DATA #####

#### plot rarefaction curves
plot(sort(sample_sums(df.prerarefaction))) #looking at sample read counts
summary(sample_sums(df.prerarefaction))

## initial plot to show rarefaction
df.matrix = as.data.frame(t(as.matrix(df.prerarefaction@otu_table)))
rarecurve(df.matrix, step=50, cex=0.5)

## sort the samples
sort(sample_sums(df.prerarefaction))

## check how many samples you'll loose at the rarefaction step
rare.threshold = which(sample_sums(df.prerarefaction) < 9862) 
rare.threshold
view(df.prerarefaction@sam_data)

## make csv of samples lost if rarefying at 9862
samples.lost <- prune_samples(sample_sums(df.prerarefaction) <= 9862, df.prerarefaction)
metadata.lost = as.data.frame(samples.lost@sam_data)
write.csv(metadata.lost, "rarefaction_threshold.csv")


## load the package stats
require(stats)
## set seed to tell R which set of random numbers to use
set.seed(5)
## check that the seed setting worked. This should return the same set of random numbers. If they don't, the seed isn't set. 
sample(10)
sample(10)

## rarefy every sample to a set number of reads (9862) here  
df.rarefied <- rarefy_even_depth(df.prerarefaction, sample.size = 9862) 



##### CHECK THAT RAREFACTION WORKED #####
## function to calculate counts along a row
abundance = function(x){
  return(sum(x,na.rm=TRUE))}

otu = as.data.frame(df.rarefied@otu_table)
## check the data are formatted correctly (see requirements for abundance function above). 
otu=as.matrix(otu)
# transpose to have sample_id as row and OTU as column (might not need to)
otu = t(otu)
otu = as.data.frame(otu)

otu$abundance = apply(otu,1,abundance)

summary(otu$abundance) 


## save rarefied data
write_rds(df.rarefied, "PAE_minocycline_rarefied_df.RDS")

