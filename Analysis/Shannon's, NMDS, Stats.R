##### WORKSPACE SETUP #####
library(phyloseq)
library(vegan)
library(stats)
library(ggplot2)
library(tidyverse)
#library(tidyr)
library(ggsignif)
library(sperich)
library(pairwiseAdonis)
library(dunn.test)
library(car)



## tell R where to get the files
setwd("C:/Users/tonym/Desktop/Year 4/Term 2/BIOL 403+503/Project")
#?import

############################VARIABLES#########################################

PAE_data = readRDS("PAE_minocycline_rarefied_df.RDS")
PAE_data
PAE_data_male <-subset_samples(PAE_data, sex=="Male")
PAE_data_female <-subset_samples(PAE_data, sex=="Female")

## GET DATA OUT

metadata = as.data.frame(as.matrix(PAE_data@sam_data))
otu = as.data.frame(t(as.matrix(PAE_data@otu_table)))
otu_table1 <- otu_table(otu, taxa_are_rows = FALSE)

# MERGE METADATA & OTU TABLE 
# the metadata and otu table together
## by can be a column but b=0 joins by rownames.
metaasv = merge(metadata, otu, by=0)

##SUBSETTING 

males_only_asv <- subset(metaasv, sex=="Male")
females_only_asv <- subset(metaasv, sex=="Female")

PAE_male <- males_only_asv[,-ncol(males_only_asv)]
PAE_meta_col <-c(1,2,3,4,5,6,7,8,9,10,11,12)
PAE_otu_col <-13:304 

################SHANNON'S INDEX ############################
metadata$obs_richness = estimate_richness(otu_table1, split = TRUE, measures = c("Observed"))
metadata$alpha <- diversity(otu, index = "shannon")


males_only = subset(metadata, sex=="Male")
females_only = subset(metadata, sex=="Female")

##########VIEW TABLES###############
#view(PAE_otu_col)
#View(otu_table1)
#View(metaasv)
#View(metadata)

############# PLOTTING SHANNON'S INDEX ALPHA DIVERSITY #####################################

#both sexes
both_sexes_graph <- ggplot(metadata, aes(x = prenatal_group, y = alpha, fill=drug_group)) + scale_fill_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  geom_boxplot(outlier.shape=NA) + xlab("Condition") +ylab("Shannon's Index") + geom_jitter(alpha = 0.4)+
  labs(title = " both females and males") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +theme(text = element_text(size = 20)) + theme(axis.text.x = element_text(size = 27))
#both_sexes_graph

#males - sex
males_only_graph <-ggplot(males_only, aes(x = prenatal_group, y = alpha, fill=drug_group)) + xlab("Condition") +ylab("Shannon's Index") + scale_fill_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  labs(title = " Males Only") + geom_boxplot(outlier.shape=NA)  + geom_jitter(alpha = 0.4)+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +theme(text = element_text(size = 20)) + theme(axis.text.x = element_text(size = 27))
#males_only_graph

#females - sex
females_only_graph <-ggplot(females_only, aes(x = prenatal_group, y = alpha, fill=drug_group)) +  scale_fill_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  xlab("Condition") +ylab("Shannon's Index")  + geom_boxplot(outlier.shape=NA)  + geom_jitter(alpha = 0.4)+
  labs(title = "Females Only") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +theme(text = element_text(size = 20)) + theme(axis.text.x = element_text(size = 27))
#females_only 

pdf("Shannon's Index (α) Plots.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "letter")          # Paper size

both_sexes_graph
males_only_graph
females_only_graph

# Closing the graphical device
dev.off() 


#ANOVA & Tukey Post-Hoc Test

a1 = aov( alpha~ prenatal_group * drug_group, data=metadata)
summary(a1)
TukeyHSD(a1)
leveneTest(alpha~ prenatal_group * drug_group, data=metadata)
plot(a1)

a2 = aov( alpha~ prenatal_group * drug_group, data=males_only)

summary(a2)
TukeyHSD(a2)
leveneTest(alpha~ prenatal_group * drug_group, data=males_only)
plot(a2)

a3 = aov( alpha~ prenatal_group * drug_group, data=females_only)
summary(a3)
TukeyHSD(a3)
leveneTest(alpha~ prenatal_group * drug_group, data=females_only)
plot(a3)


#################
# MDS Functions #
#################
# Function to run the ellipse calculation for each of column specified by group
# Uasge: calculate_ellipse(seagrass,"region")
# Will output a data frame with three columns,
#    one called "filter" for the column you are filtering by
#    and one each for MDS1 and MDS2
calculate_ellipse = function(df,filter_by){
  # Function to calculate the 95% ci ellipses for each sample
  veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ellipse <- data.frame()
  for(g in unique(df[[filter_by]])){
    df_ellipse <- rbind(df_ellipse,cbind(as.data.frame(with(df [df[filter_by]==g,],
                                                            veganCovEllipse(cov.wt(cbind(MDS1,MDS2),
                                                                                   wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2))))),filter=g))
  }
  return(df_ellipse)
}

# Function to test stress levels at different
# Usage: checkMDSdim(seagrass)
# Optional arguments include
#    iter: the number of times to run MDS for each dimension
#    dimension: a vector of dimensions to run
#    distance, try, and trymax as inputs of metaMDS
checkMDSdim = function(data,iter=1,dimensions=1:5,distance="bray",try=20,trymax=20){
  n = length(dimensions*iter)
  stress = rep(0,n)
  c = 0
  for (d in 1:n){
    temp = rep(0,iter)
    for (i in 1:iter){
      c = c + 1
      stress[c] = metaMDS(data,k=d,distance=distance,try=try,trymax=trymax)$stress
    }
  }
  output = data.frame(dimension = rep(dimensions,each=iter),
                      stress = stress)
  return(output)
}

combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

################################BETA DIVERSITY - BRAY-CURTIS ###########################################
highlight_var  = "drug_group"
variable_to_filter_by2 = "prenatal_group"

# BOTH SEXES --------------------------------------------------------------------------------------
#{r, Plotting Stress values across Dimensions, include=TRUE, echo=FALSE}
PAE_MDSdim_both<-checkMDSdim((metaasv[-PAE_meta_col]),iter = 5)
plot_stress_both = ggplot(PAE_MDSdim_both, aes(dimension, stress))+
  geom_point(size = 2)+
  theme_classic()+
  theme(text = element_text(size=15,face = "bold"))

#{r, NMDS plot, include=TRUE, echo=FALSE }
PAE_MDS_both <- metaMDS(metaasv[-PAE_meta_col])
PAE_MDSmeta_both = data.frame(cbind(metaasv[,PAE_meta_col],PAE_MDS_both$points))
plot_MDS_both = ggplot(PAE_MDSmeta_both, aes(x = MDS1,y = MDS2, color=drug_group, shape=prenatal_group))+
  geom_point(aes_string(shape=variable_to_filter_by2, color=highlight_var), size= 5) + scale_colour_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  coord_fixed()+ labs(title = " Figure 3 : Both ")+
  theme_classic()+
  theme(text = element_text(size = 15, face = "bold")) 
plot_MDS_both 

PAE_MDS_both

# MALES --------------------------------------------------------------------------------------
#{r, Plotting Stress values across Dimensions, include=TRUE, echo=FALSE}
PAE_MDSdim_males<-checkMDSdim((males_only_asv[-PAE_meta_col]),iter = 5)
plot_stress = ggplot(PAE_MDSdim_males, aes(dimension, stress), include=TRUE, echo=FALSE)+
  geom_point(size = 2)+ 
  theme_classic()+
  theme(text = element_text(size=15,face = "bold"))

#{r, NMDS plot, include=TRUE, echo=FALSE }
PAE_MDS_males <- metaMDS(males_only_asv[-PAE_meta_col])
PAE_MDSmeta_males = data.frame(cbind(males_only_asv[,PAE_meta_col],PAE_MDS_males$points))
plot_MDS_males = ggplot(PAE_MDSmeta_males, aes(x = MDS1,y = MDS2))+
  geom_point(aes_string(shape=variable_to_filter_by2, color=highlight_var), size= 5) + scale_colour_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  coord_fixed()+ labs(title = " Figure 1: Males Only") +
  theme_classic()+
  theme(text = element_text(size = 15, face = "bold"))
plot_MDS_males


# FEMALES --------------------------------------------------------------------------------------
#{r, Plotting Stress values across Dimensions, include=TRUE, echo=FALSE}
PAE_MDSdim_females<-checkMDSdim((females_only_asv[-PAE_meta_col]),iter = 5)
plot_stress_females = ggplot(PAE_MDSdim_females, aes(dimension, stress))+
  geom_point(size = 2)+ 
  theme_classic()+
  theme(text = element_text(size=15,face = "bold"))

#{r, NMDS plot, include=TRUE, echo=FALSE }
PAE_MDS_females <- metaMDS(females_only_asv[-PAE_meta_col])
PAE_MDSmeta_females = data.frame(cbind(females_only_asv[,PAE_meta_col],PAE_MDS_females$points))
plot_MDS_females = ggplot(PAE_MDSmeta_females, aes(x = MDS1,y = MDS2))+
  geom_point(aes_string(shape=variable_to_filter_by2, color=highlight_var), size= 5)+ scale_colour_manual(values = c("#b220d1", "#e87e00", "#31FFE8"))+
  coord_fixed()+ labs(title = " Figure 2: Females Only") +
  theme_classic()+
  theme(text = element_text(size = 15, face = "bold"))
plot_MDS_females

pdf("MDS Plots Bray-Curtis.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "letter")          # Paper size

plot_MDS_both
plot_MDS_males
plot_MDS_females


# Closing the graphical device
dev.off() 

######################### SUBSETTING METAASV ##################
males_only_asv <- subset(metaasv, sex=="Male")
females_only_asv <- subset(metaasv, sex=="Female")

##Run the PERMANOVA
p1 = adonis2(metaasv[,-c(1:12)]~prenatal_group * drug_group, data=metaasv)
print(p1)

p2 = adonis2(males_only_asv[,-c(1:12)]~prenatal_group * drug_group, data=males_only_asv)
print(p2)

p3 = adonis2(females_only_asv[,-c(1:12)]~prenatal_group * drug_group, data=females_only_asv)
print(p3)

#########################BETA-DISPERSION TESTS####################################

project_bray.rarefied <- phyloseq::distance(PAE_data, method = "bray")
project_bray.rarefied.male <- phyloseq::distance(PAE_data_male, method = "bray")
project_bray.rarefied.female <- phyloseq::distance(PAE_data_female, method = "bray")

sample_df <- data.frame(sample_data(metaasv))
sample_df1 <- data.frame(sample_data(males_only_asv))
sample_df2 <- data.frame(sample_data(females_only_asv))

beta.both_prenatal_group <- betadisper(project_bray.rarefied, sample_df$prenatal_group) 
b1.1 = permutest(beta.both_prenatal_group) 
b1.1
beta.both_drug_group <- betadisper(project_bray.rarefied, sample_df$drug_group) 
b1.2 = permutest(beta.both_drug_group) 
b1.2

beta.male_prenatal_group <- betadisper(project_bray.rarefied.male, sample_df1$prenatal_group) 
b1.3 = permutest(beta.male_prenatal_group) 
b1.3
beta.male_drug_group <- betadisper(project_bray.rarefied.male, sample_df1$drug_group) 
b1.4 = permutest(beta.male_drug_group) 
b1.4

beta.female_prenatal_group <- betadisper(project_bray.rarefied.female, sample_df2$prenatal_group) 
b1.5 = permutest(beta.female_prenatal_group) 
b1.5
beta.female_drug_group <- betadisper(project_bray.rarefied.female, sample_df2$drug_group) 
b1.6 = permutest(beta.female_drug_group) 
b1.6

################################### PERMANOVA POST-HOC #########################################

metasv_control<- subset(metaasv, prenatal_group=="Control")
metasv_ethanol<- subset(metaasv, prenatal_group=="Ethanol")

males_only_control<- subset(metaasv, sex=="Male" & prenatal_group=="Control")
males_only_ethanol <- subset(metaasv, sex=="Male" & prenatal_group=="Ethanol")
females_only_control <- subset(metaasv, sex=="Female" & prenatal_group=="Control")
females_only_ethanol <- subset(metaasv, sex=="Female"& prenatal_group=="Ethanol")

## Conduct a post-hoc test

pair.id = pairwise.adonis(metasv_control[,-c(1:12)], factors = (metasv_control$drug_group))
print(pair.id)

pair.id = pairwise.adonis(metasv_ethanol[,-c(1:12)], factors = (metasv_ethanol$drug_group) )
print(pair.id)

pair.id = pairwise.adonis(males_only_control[,-c(1:12)], factors = (males_only_control$drug_group))
print(pair.id)

pair.id = pairwise.adonis(males_only_ethanol[,-c(1:12)], factors = (males_only_ethanol$drug_group) )
print(pair.id)

pair.id = pairwise.adonis(females_only_control[,-c(1:12)], factors = (females_only_control$drug_group) )
print(pair.id)

pair.id = pairwise.adonis(females_only_ethanol[,-c(1:12)], factors = (females_only_ethanol$drug_group) )
print(pair.id)




