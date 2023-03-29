library(tidyverse)
library(phyloseq)
library(indicspecies)
library(DESeq2)

setwd("~/Desktop/PAE")

## read in files
df = readRDS("PAE_minocycline_filtered_phyloseq.RDS")
df
df_minocycline <- subset_samples(df, drug_group == "No drug")

## First comparing control group without minocycline and PAE group without minocycline
minocycline_DESeq <- phyloseq_to_deseq2(df_minocycline, design = ~prenatal_group)
minocycline_DESeq <- DESeq(minocycline_DESeq, test = "Wald", fitType = "parametric")

minocycline_DESeq_results <- results(minocycline_DESeq)
minocycline_DESeq_results <- data.frame(cbind(minocycline_DESeq_results, df_minocycline@tax_table))

plot <- ggplot(minocycline_DESeq_results, aes(x = Rank6, y = log2FoldChange, color = Rank6))+
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
plot

## Filter for signif p values p<0.05
minocycline_DESeq_results_sig <- filter(minocycline_DESeq_results,padj<0.05)
minocycline_DESeq_results_narm <- filter(minocycline_DESeq_results)

# function to combine p-values
combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

# mean of each phyla
minocycline_DEseq_results_phyla_summary<-minocycline_DESeq_results_sig %>%
  group_by(Rank2) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
            padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue))

# mean of each genus
minocycline_DEseq_results_genus_summary<-minocycline_DESeq_results_sig %>%
  group_by(Rank6) %>%
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
genus_plot <- ggplot(minocycline_DEseq_results_genus_summary,aes(x=reorder(Rank6, -mean_2_fold_change),y=mean_2_fold_change,color = Rank6))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "none",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
  labs(x="Genus",y="Log2 Fold Change",color= "Genus")


## Subset df comparing control group without minocycline and PAE group with minocycline
df_ctrl <- subset_samples(df, prenatal_group == "Control" & drug_group == "No drug")
df_PAE <- subset_samples(df, prenatal_group == "Ethanol" & drug_group != "No drug")
df_hyp2 <- merge_phyloseq(df_ctrl, df_PAE)

hyp2_DESeq <- phyloseq_to_deseq2(df_hyp2, design = ~prenatal_group)
hyp2_DESeq <- DESeq(hyp2_DESeq, test = "Wald", fitType = "parametric")

hyp2_DESeq_results <- results(hyp2_DESeq)
hyp2_DESeq_results <- data.frame(cbind(hyp2_DESeq_results, df_hyp2@tax_table))

hyp2_DESeq_results_sig <- filter(hyp2_DESeq_results,padj<0.05)
hyp2_DESeq_results_narm <- filter(hyp2_DESeq_results)

hyp2_DEseq_results_phyla_summary<-hyp2_DESeq_results_sig %>%
  group_by(Rank2) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
                   padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue))

hyp2_DEseq_results_genus_summary<-hyp2_DESeq_results_sig %>%
  group_by(Rank6) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange, na.rm = TRUE),
                   padj=combine_pvalues(padj),countmean = mean(baseMean, na.rm = TRUE),pvalue=combine_pvalues(pvalue)) %>%
  filter(Rank6 != "uncultured")

phyla2_plot <- ggplot(hyp2_DEseq_results_phyla_summary,aes(x=reorder(Rank2, -mean_2_fold_change),y=mean_2_fold_change,color = Rank2))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "none",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
  labs(x="Phyla",y="Log2 Fold Change",color= "Phyla")

genus2_plot <- ggplot(hyp2_DEseq_results_genus_summary,aes(x=reorder(Rank6, -mean_2_fold_change),y=mean_2_fold_change,color = Rank6))+
  geom_point(size=3.5) +
  geom_hline(yintercept = 0,size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold",size=12),
        legend.position = "none",
        axis.text.x.bottom  = element_text(size=9,angle = -45)) +
  labs(x="Genus",y="Log2 Fold Change",color= "Genus")
