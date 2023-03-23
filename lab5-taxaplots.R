##### SET UP #####
## load packages
library(phyloseq)
library(tidyverse)
library(plyr)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text.y = element_text(size = 12, colour = "black"),
                                    axis.title = element_text(size=15, face="bold"),
                                    strip.text = element_text(color="black", size=10),
                                    legend.text=element_text(size=10),
                                    axis.line = element_line(colour = "black"),
                                    axis.text.x = element_blank(),))

## load function
dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  ## if your metadta is empty after running this, you need to use otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_id")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## set working directory
setwd("~/Desktop/PAE")
filepath ="~/Desktop/PAE"

## read in the data
PAE.phyloseq = read_rds("PAE_minocycline_filtered_phyloseq.RDS")

## get the number of reads in each sample in your phyloseq object
PAE.phyloseq@sam_data$number_of_reads = sample_sums(PAE.phyloseq)

## subset the phyloseq object to only keep the samples you want
PAE.phyloseq = subset_samples(PAE.phyloseq, prenatal_control = "Control")

## get the PAE data out of phyloseq
## I always call this dataframe mot becayse it's the Metadata, OTU table, and Taxonomy, in that order
mot = dephyloseq(PAE.phyloseq)

## get realtive abundance of otu
mot$relative_abundance = as.numeric(mot$asv_abundance)/as.numeric(mot$number_of_reads)

## make taxa names for plot column
mot$plot_names = paste0(mot$Rank4, "; ", mot$Rank6)

##### LOOP TO SELECT TOP 15 GENERA IN SAMPLES ######
## make variable to track sample types and salinities
mot$group = paste0(mot$prenatal_group, "-", mot$drug_group)


## summarize mto
mot.sum = ddply(mot, c("group", "plot_names"),
                summarise,
                N=length(relative_abundance),
                sum = sum(relative_abundance))


## get list of all samples
samplegroups = unique(mot.sum$drug_group)

## sort data by relative abundance
sorted = mot.sum[order(-mot.sum$sum),]

## make empty dataframe to store taxa
top.df = NULL


## start loop
for(i in samplegroups) {
  for(j in i) {
    ## subset dataframe by samples
    sample = subset(sorted, sorted$group %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:15),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}

## add identifier for top and bottom taxa
top.df$place = "top_15"

##### SET UP TAXA PLOT #####
top.df <- data.frame(top.df)
## join the top taxa and existing mto dataframe
mot.top = full_join(mot, top.df)

## make the empty "place" cells say bottom
mot.top$place = replace(mot.top$place, is.na(mot.top$place), "bottom")

## replace plot_names for bottom taxa with Other
mot.top[mot.top$place == "bottom",]$plot_names <- "Others"



##### MAKE COLOR LIST #####
## how many do we need?
taxa = unique(mot.top$plot_names)

## write a csv with the names
#write.csv(taxa, "colors_taxaplot_biol403503.csv")
colors = read.csv("colors_taxaplot_biol403503.csv")

##### MAKE TAXAPLOT #####
# Assign taxa names to colors
scolors <- colors$plot_colors
names(scolors) <- colors$plot_names

## order mot.top by decreasing relative abundance
mot.top = mot.top[order(-mot.top$relative_abundance),]

## get list of factors in order
natural.genus.order = as.list(c(unique(mot.top$plot_names)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
mot.top$plot_names = factor(mot.top$plot_names, levels=c(plot.order))

## get list to cycle through
mot.top$taxaplotgroup = paste0(mot.top$Study, "-", mot.top$Group)
taxaplot.groups = unique(mot.top$taxaplotgroup)


## field taxaplot loop
for (i in taxaplot.groups){
  for (j in i){
    sub.df = subset(mot.top, mot.top$taxaplotgroup== c(j))
    
    myplot=ggplot(sub.df, aes(x=as.character(Row.names), y=as.numeric(relative_abundance), 
                              fill=as.factor(plot_names)))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values=scolors)+
      facet_grid(.~Study+Group, scales="free", space="free")+
      labs(x=" ", y="Genus relative abundance in samples", fill="Order; Genus")+
      guides(fill=guide_legend(ncol=1))
    
    myplot
    
    ## save plot
    ggsave(myplot, filename=paste0(j,"_practice_taxaplot",".png",sep=""), width=12, height=7,
    path=paste0(filepath, "/practice_taxaplots/"))
  }
}













