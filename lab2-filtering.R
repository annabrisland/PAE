##### WORKSPACE SETUP #####
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
df = subset_taxa(df, Kingdom="Archaea"|"Bacteria")

##### FILTERING - REMOVE SAMPLES WITH LESS THAN N READS #####
# my N is usually >1000
## Remove samples with less than N reads 
sample_sums(df)
plot(sort(sample_sums(df)))

## add number of total sequences in a sample (Read depth) to the metadata 
df@sam_data$read_depth_noofftargets = sample_sums(df) 

## check which samples have less than N reads
which(df@sam_data$read_depth_noofftargets < 1000)


## remove the samples by your threshold 
df.pruned <- prune_samples(sample_sums(df) >= 1000, df)

## write file to know which samples were lost here. This is important for the methods section. 
df.below1000 <- prune_samples(sample_sums(df) < 1000, df)
df.below1000 = as.matrix(df.below1000@sam_data)
write_rds(df.below1000, "SF_samples_less_than_N.csv")

##### FILTERING - REMOVE INDIVIDUAL ASVS WITH LESS THAN N READS #####
## N is probably around 100 
## extract OTU dataframe from phyloseq object
otu.pruned <- as.data.frame(as.matrix(otu_table(df.pruned)))

## remove ASVs (rows) with less than N reads accross whole dataset but keep all samples
#!# make sure asv sequence is rownames and sample id is column name
otu.pruned$rowsum = rowSums(otu.pruned)

## remove low frequency asvs
otu.pruned = subset(otu.pruned, otu.pruned$rowsum>100)

## remove rowsum column from your OTU table 
otu = subset(otu.pruned, select=-c(rowsum))


##### FILTERING - REMOVE ASVs FOUND IN N SAMPLES OR LESS #####
## your N is probably between 2 and 5

# has sample ID as column name and asv is as row name. Needs to be this way to use richness function 
# function to calculate richness, sums along a row (OTU)
richness = function(x){return(sum(x>0))}

## calculate richness on entire dataframe
otu$richness = apply(otu,1,richness) # use all columns of otu dataframe
summary(otu$richness)

## remove OTU (rows) with richness N or less (found in two samples or less) but keep all samples (columns)
otu = subset(otu, otu$richness> 2)
## check that it worked
summary(otu$richness)
## remove richness column
otu = subset(otu, select=-c(richness))


##### DENOISING - MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES N OR LESS 0 #####
## your N is probably between 2 and 10

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

## remove samples with very highest sample counts (might not be required for your dataset)
histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)
#df.prerarefaction = subset_samples(df.prerarefaction, read_depth_filtered<75000)
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
# N here will be above your "FILTERING - REMOVE SAMPLES WITH LESS THAN N READS" N number
rare.threshold = which(sample_sums(df.prerarefaction) < 9862) 
rare.threshold
view(df.prerarefaction@sam_data)

## make csv of samples lost if rarefying at N
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

## rarefy every sample to a set number of reads (N) here  
df.rarefied <- rarefy_even_depth(df.prerarefaction, sample.size = 9862) 



##### CHECK THAT RAREFACTION WORKED #####
## function to calculate counts along a row
#!# this function needs sample_id as a row and asvs as columns
abundance = function(x){
  return(sum(x,na.rm=TRUE))}

otu = as.data.frame(df.rarefied@otu_table)
## check the data are formatted correctly (see requiermtns for abundance function above). 
otu=as.matrix(otu)
# transpose to have sample_id as row and OTU as column (might not need to)
otu = t(otu)
otu = as.data.frame(otu)

otu$abundance = apply(otu,1,abundance)
## should be your rarefaction number 9N) exactly
summary(otu$abundance) 


## save rarefied data
write_rds(df.rarefied, "PAE_minocycline_rarefied_df.RDS")

