
# already installed packages:
library(UpSetR)
library(plyr)
library(limma)
library(UpSetR)
library(reshape2)
library(dplyr)
library(tidyr)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(edgeR)
library(ggupset)
library(DOSE)
library(magrittr)
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

BiocManager::install("org.Hs.eg.db")
BiocManager::install("GenomeInfoDb")

library(DESeq2)

library(GenomeInfoDb)
install.packages("xtable")


## start code ----
all_featureCounts <- read.table(file = "FeatureCounts/complete_sample_FC_final_analysis_0.2.txt", 
                                sep = "\t", 
                                stringsAsFactors = FALSE,
                                header = TRUE, 
                                row.names = 1)
all_featureCounts <- all_featureCounts[ -c(1:5) ] # remove first 5 columns
names(all_featureCounts)
colnames(all_featureCounts)[1:21] <-c("primary_m1e2_4P", "adrenal_m1e2_22", "bone_m1e2_2_II", 
                                      "primary_m2e2_5P", "adrenal_m2e2_12_I", "brain_m2e2_E", "bone_m2e2_5_I",
                                      "primary_m3e2_1P", "adrenal_m3e2_2", "kidney_m3e2_12", "tibia_m3e2_C", "femur_m3e2_15_I",
                                      "primary_m4ed_2P", "bone_m4ed_16",
                                      "primary_m5ed_3P", "brain_m5ed_11_II",
                                      "primary_m6ed_25", "bone_m6ed_10", "brain_m6ed_6",
                                      "primary_m7ed_26", "brain_m7ed_8_I")
names(all_featureCounts)

## subset ED bone and primary samples (n=2) apply BCV from E2 bone samples (n=3) of 0.3
# see if similar increase in ESR1 expresion as found in the ED brain sample group

counts <- all_featureCounts[c(20,21)] 
names(counts)
Groups <- c(1,2) # 1=primary 2=mets

# DGEList object ----
DGE <- DGEList(counts = counts, group = Groups, genes = row.names(counts))
dim(DGE)

# Filter out low expression genes (fewer than  CPM) in two or more samples ----
DGE <- DGE[rowSums(cpm(DGE) > 1) >= 2, keep.lib.sizes = FALSE] 
dim(DGE) 

# TMM Normalize ----
DGE <- calcNormFactors(DGE)
DGE$samples

# lib size bar plot ----
barplot(DGE$samples$lib.size,names=colnames(DGE),las=2)
#title("2P.25v16.10")

# estimate dispersions for singe paired samples---
et <- exactTest(DGE, dispersion = 0.3) # katia code  
topTags(et, n=50) 
results <- as.data.frame(topTags(et, n = dim(DGE)[1]))
#write.csv(results_25_10, file = "1P primary v 15_I andc E2/RESULTS_E2_1P_15.1.csv", row.names = FALSE)

# # up/down DEGs ----
sig <- decideTestsDGE(et, p = 0.05, adjust = "none")
summary(sig)

# filter sig results ----
#results_25_10 <- as.data.frame(topTags(et_ed_bone, n = dim(DGE_ed_bone)[1]))
sig_results <- results[abs(results$logFC) > 1 & results$PValue < 0.05,  ]
dim(sig_results)
write.csv(sig_results, file = "EdgeR_output_files/SIG_RESULTS_ED_ZR751_26v8_1.csv", row.names = FALSE)


install.packages("DESeq2")

#transform count data
rld<-rlog(deg)





# Volcano Plot ----
plot(-log10(results_25_10$PValue) ~ results_25_10$logFC, 
     xlab = "log2 FC", ylab = "-log10 p-value", 
     main = "4P MCF-7 E2 primary v matched 2.II tibia")
points(-log10(results_25_10$PValue[row.names(results_25_10) %in% sig_ED_25_10$genes]) ~ 
         results_25_10$logFC[row.names(results_25_10) %in% sig_ED_25_10$genes], col = "red", pch = 16)

## GO & KEGG pathway analysis ----
library(GO.db)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

GO <- goana(et_ed_bone, species="Hs")
topGO(GO, sort = "up")
write.csv(GO, file = "25 primary V 10 femur met/GO_ED_25_10.csv", row.names = FALSE)

kegg <- kegga(et_ed_bone, species="Hs")
topKEGG(kegg, sort = "up")
write.csv(kegg, file = "25 primary V 10 femur met/KEGG_ED_25_10.csv", row.names = FALSE)

## WIKI GSEA ----


wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# upload data and format: data(geneList) # data set structure reference ----
d <- read.csv(file = "upsetplot input/filter out lowly expressed by 2/gene.list.csv") # col 1 = entrez ID, col 2 = logFC 
mygenelist <- d[,2] ## feature 1: numeric vector
names(mygenelist) <- as.character(d[,1]) ## feature 2: named vector
mygenelist <- sort(mygenelist, decreasing = TRUE) ## feature 3: decreasing order
de <- names(mygenelist)[abs(mygenelist) > 1] # subset data to genes with FC greater than 1
head(de)
#ewp <- enricher(de, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
#head(ewp)
#write.csv(ewp, file = "2P primary V 16 humerus met/ew0WIKI_GSEA_READABLE_ED_2P_16.csv", row.names = FALSE)

ewp2 <- GSEA(mygenelist, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff=1)
head(ewp2)
write.csv(ewp2, file = "upsetplot input/filter out lowly expressed by 2/WIKI_GSEA_READABLE_E2_4P_2.II.csv", row.names = FALSE)

#alternative code; use function gseNCH ----
#edo2 <- gseNCG(mygenelist, nPerm=10000)
#p8 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA") #******PLOTTTTTT
#plot(p8)

p4 <- dotplot(ewp2, showCategory=30) + ggtitle("WIKI GSEA: (E2) MCF-7 4P primary v. 2_II tibia")
plot(p4)
#plot_grid(p1, p2, ncol=2)

## convert gene ID to Symbol
edox <- setReadable(ewp2, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=mygenelist)
plot(p1)


## categorySize can be scaled by 'pvalue' or 'geneNum'

p2 <- cnetplot(edox, categorySize="pvalue", foldChange=mygenelist) 
plot(p2)
p3 <- cnetplot(edox, foldChange=mygenelist, circular = TRUE, colorEdge = TRUE)#### should contain at least two columns!!!
plot(p3)
p3
#cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p7 <- heatplot(ewp2, foldChange = mygenelist)
plot(p7)
p7 + theme(axis.text.x = element_text(angle = 90))



#   -------- FC subset code, if no sig data ------

# filter log FC >2 results ----
results_25_10_FC<- results_25_10[(results_25_10$logFC)> 1,  ]
dim(results_25_10_FC)
write.csv(results_25_10_FC, file = "25 primary V 10 femur met/FC_results_25v10_FC_gt_1.csv", row.names = FALSE)




