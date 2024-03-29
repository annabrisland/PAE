library(tidyverse)
library(phyloseq)
library(indicspecies)
library(DESeq2)

setwd("~/Desktop/PAE")

# function to combine p-values
combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

## read in files
df = readRDS("PAE_minocycline_filtered_phyloseq.RDS")

##Agglomerate taxa of the same type at phylum level
df_glom_phy <- tax_glom(df, taxrank = "Rank2")

##Agglomerate taxa of the same type at genus level
df_glom <- tax_glom(df, taxrank = "Rank6")

df_minocycline <- subset_samples(df_glom_phy, drug_group == "No drug")

## First comparing control group without minocycline and PAE group without minocycline
# Create DESeq object and run DESeq
minocycline_DESeq <- phyloseq_to_deseq2(df_minocycline, design = ~prenatal_group)
minocycline_DESeq <- DESeq(minocycline_DESeq, test = "Wald", fitType = "parametric")

minocycline_DESeq_results <- results(minocycline_DESeq)
minocycline_DESeq_results <- data.frame(cbind(minocycline_DESeq_results, df_minocycline@tax_table))

# Plot all results
plot <- ggplot(minocycline_DESeq_results, aes(x = Rank6, y = log2FoldChange, color = Rank6))+
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))

## Filter for signif p values p<0.05
minocycline_DESeq_results_sig <- filter(minocycline_DESeq_results,padj<0.05)

# mean of each phyla
minocycline_DEseq_results_phyla_summary<-minocycline_DESeq_results_sig %>%
  group_by(Rank2) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
            padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue))

# mean of each genus
minocycline_DEseq_results_genus_summary<-minocycline_DESeq_results_sig %>%
  group_by(Rank6, Rank2) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
                   padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue)) %>%
  filter(Rank6 != "uncultured")


## Plot phyla
phyla_plot <- ggplot(minocycline_DEseq_results_phyla_summary,aes(x=reorder(Rank2, -mean_2_fold_change),y=mean_2_fold_change,color = Rank2))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "none",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
        labs(x="Phylum",y="Log2 Fold Change",color= "Phylum")

## Plot genus
genus_plot <- ggplot(minocycline_DEseq_results_genus_summary,aes(x=reorder(Rank6, -mean_2_fold_change),y=mean_2_fold_change,color = Rank2))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "right",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
  labs(x="Genus",y="Log2 Fold Change",color= "Phyla")

### Comparing PAE with minocycline vs PAE without minocycline
df_ethanol <- subset_samples(df_glom, prenatal_group == "Ethanol")

# Create DESeq object and run DESeq
ethanol_DESeq <- phyloseq_to_deseq2(df_ethanol, design = ~drug_group)
ethanol_DESeq <- DESeq(ethanol_DESeq, test = "Wald", fitType = "parametric")

ethanol_DESeq_results <- results(ethanol_DESeq)
ethanol_DESeq_results <- data.frame(cbind(ethanol_DESeq_results, df_ethanol@tax_table))

ethanol_DESeq_results_sig <- filter(ethanol_DESeq_results,padj<0.05)
ethanol_DESeq_results_narm <- filter(ethanol_DESeq_results)

ethanol_DEseq_results_genus_summary<-ethanol_DESeq_results_sig %>%
  group_by(Rank6, Rank2) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
                   padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue)) %>%
  filter(!Rank6 %in% c("uncultured", "uncultured_bacterium"))

# Plot significant genera
genus_ethanol_plot <- ggplot(ethanol_DEseq_results_genus_summary,aes(x=reorder(Rank6, -mean_2_fold_change),y=mean_2_fold_change,color = Rank2))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  scale_color_manual(values=c("#FF69B4", "#FF8C00", "#4682B4")) +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "right",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
  labs(x="Genus",y="Log2 Fold Change",color= "Phyla")


