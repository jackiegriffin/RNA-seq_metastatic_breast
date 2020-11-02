BiocManager::install('DOSE')
BiocManager::install('clusterProfiler')
BiocManager::install('tibble')

if (!require("DOSE")) install.packages("DOSE")
library(DOSE)
if(!require("clusterProfiler")) install.packages("clusterProfiler")
library(clusterProfiler)
if (!require("enrichplot")) install.packages("enrichplot")
library(enrichplot)
if (!require("forcats")) install.packages('forcats')
library('forcats')
if (!require("ggplot2")) install.packages("ggplot2")
library('ggplot2')

if (!require("ggstance")) install.packages("ggstance")
library('ggstance')
if (!require("plyr")) install.packages("plyr")
library('plyr')
if (!require("dplyr")) install.packages("dplyr")
library('dplyr')
if (!require("viridis")) install.packages("viridis")
library("viridis")
devtools::install_github("tidyverse/dplyr")
library(rlang)
library(tibble)
library(tidyverse)

NES_filtered_B<-read.csv("25 primary V 10 femur met/split_plotWIKI_GSEA_READABLE_ED_25_10_FC_gt1.csv", header = TRUE, stringsAsFactors = FALSE)
#split_NES_filtered_B <- arrange(NES_filtered_B, abs(nes)) %>% 
 # group_by(sign(nes)) %>% 
  #slice(1:30)

ggplot(NES_filtered_B, aes(nes, fct_reorder(Description, nes), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL) + ggtitle("Pathway dysregulation in ZR751 femur metastasis (ED)", 
                                         subtitle = "Vs. matched primary (ED)") +
 xlab("Normalized Enrichment Score") + labs(caption = "(ED) ZR75-1 25 primary (ERa-) & 10 femur (ERa+)") +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  theme(axis.title = element_text(size = 12)) +
  theme(title = element_text(size = 14)) +
  theme(plot.subtitle = element_text(size = 12)) +
  theme(plot.caption = element_text(size = 10))



# --------------------------------------------------------------------------
#net10<-read.csv("networkin_input.csv", header = TRUE, stringsAsFactors = FALSE)
ggplot(NES_filtered_B, showCategory = 50, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=nes, size = count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=FALSE), direction = -1) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("gene ratio (data/signature set)") +
  ylab(NULL) + 
  labs(caption = " ") +
  theme(axis.text = element_text(size = 12, colour = 'black')) +
  theme(axis.title = element_text(size = 12)) +
  theme(title = element_text(size = 14)) +
  theme(plot.subtitle = element_text(size = 12)) +
  ggtitle("Pathway dysregulation in ZR751 femur metastasis (ED)", subtitle = "Vs. matched primary (ED)") +
  theme(plot.caption = element_text(size = 10))



