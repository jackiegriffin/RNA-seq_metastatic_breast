# libraries ----
  library(plyr)
  library(UpSetR)
  library(reshape2)
  library(dplyr)
  library(tidyr)
  library(gplots)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(edgeR)
  library(limma)
  library(ggupset)
  library(DOSE)
  library(magrittr)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(clusterProfiler)
  library(GenomeInfoDb)
  library(xtable)

# load featurecounts & metadata ----
  counts <- read.table(file = "2.12.2021/EdgeR/Data_input/feature_counts.txt", 
                                  sep = "\t", 
                                  stringsAsFactors = FALSE,
                                  header = TRUE, 
                                  row.names = 1)
  metadata <- read.csv(file = "2.12.2021/EdgeR/Data_input/metadata_complete.csv", 
                                  stringsAsFactors = FALSE, header = FALSE)

# format df ----
  colnames(counts)
  head(metadata)
  counts_filt <- counts[ -c(1:5) ]
  colnames(counts_filt)
  
  colnames(counts_filt)[1:21] <-c("primary_m1e2_4P", "adrenal_m1e2_22", "bone_m1e2_2_II", 
                                        "primary_m2e2_5P", "adrenal_m2e2_12_I", "brain_m2e2_E", "bone_m2e2_5_I",
                                        "primary_m3e2_1P", "adrenal_m3e2_2", "kidney_m3e2_12", "tibia_m3e2_C", "femur_m3e2_15_I",
                                        "primary_m4ed_2P", "bone_m4ed_16",
                                        "primary_m5ed_3P", "brain_m5ed_11_II",
                                        "primary_m6ed_25", "bone_m6ed_10", "brain_m6ed_6",
                                        "primary_m7ed_26", "brain_m7ed_8_I")
# COMPLETE DATASET ----
# Asign groups
  Groups_complete <- c(1,2,2, # 1 = primary 2 = met
              1,2,2,2,
              1,2,2,2,2,
              1,2,
              1,2,
              1,2,2,
              1,2) 
# DGEList
  DGE_complete <- DGEList(counts = counts_filt, group = Groups_complete, genes = row.names(counts_filt))
  dim(DGE_complete)
  
# Filter
  DGE_complete_filt <- DGE_complete[rowSums(cpm(DGE_complete) > 0.1) >= 2, keep.lib.sizes = FALSE] 
  dim(DGE_complete_filt)
  
# Normalize (TMM)
  DGE_complete_filt <- calcNormFactors(DGE_complete_filt, method = "TMM")
  DGE_complete_filt$samples

  
# counts & MDS plots 
  cols <- as.numeric(DGE_complete_filt$samples$group)+2
  par(mfrow=c(2,1))
  barplot(colSums(DGE_complete_filt$counts), las=2, main="Counts per index",
            col=cols, cex.names=0.5, cex.axis=0.8)
  legend("topright", legend=c("Primary", "Metastases"), col=c(3,4), pch=15)
  
  cols2 <- DGE_complete_filt$samples$Groups_complete
  plotMDS(DGE_complete_filt, col=cols, main="MDS Plot")
  legend("bottomright", legend=c("Primary", "metastases"), col=c(3,4), pch=15)
  
  

  
  
# ED SUBSET to estimate dispersion for single rep data ----
     ED_counts <- counts_filt[c(15:17, 19:21)] 
     colnames(ED_counts)
   
  # Asign groups
    Groups_ED_subset <- c(1,2,1,2,1,2) # 1 = primary 2 = met
                       
  # DGEList
    DGE_ED_subset <- DGEList(counts = ED_counts, group = Groups_ED_subset, genes = row.names(ED_counts))
    dim(DGE_ED_subset)
    
  # Filter
    DGE_ED_subset <- DGE_ED_subset[rowSums(cpm(DGE_ED_subset) > 0.1) >= 6, keep.lib.sizes = FALSE] 
    dim(DGE_ED_subset)
  
  # Normalize (TMM)
    DGE_ED_subset <- calcNormFactors(DGE_ED_subset, method = "TMM")
    DGE_ED_subset$samples
  
       # counts & MDS plots ----
        # cols <- as.numeric(DGE_ED_subset$samples$group)+2
        # par(mfrow=c(2,1))
        # barplot(colSums(DGE_ED_subset$counts), las=2, main="Counts per index",
        #         col=cols, cex.names=0.5, cex.axis=0.8)
        # legend("topright", legend=c("Primary", "Metastases"), col=c(3,4), pch=15)
        # 
        # cols2 <- DGE_ED_subset$samples$Groups_ED_subset
        # plotMDS(DGE_ED_subset, col=cols, main="MDS Plot")
        # legend("bottomleft", legend=c("Primary", "metastases"), col=c(3,4), pch=15)
    
    

# BCV ----
    DGE_ED_subset <- estimateCommonDisp(DGE_ED_subset, verbose = TRUE) # BCV = 0.47
    
   
# sINGLE REPLICATE DF's ----
  # Assign groups 
    counts <- all_featureCounts[c(20,21)] 
    names(counts)
    Groups <- c(1,2) # 1 = primary 2 = met

  # DGEList
    DGE <- DGEList(counts = counts, group = Groups, genes = row.names(counts))
    dim(DGE)
    
 # Filter 
    DGE <- DGE[rowSums(cpm(DGE) > .1) >= 2, keep.lib.sizes = FALSE] 
    dim(DGE)
    
  # Normalize (TMM)
    DGE <- calcNormFactors(DGE, method = "TMM")
    DGE$samples
    plotMD(cpm(DGE, log = TRUE), column = 1) # check performance
    abline(h=0, col="red", lty=2, lwd=2)
    
  # glmFit ---- (BCV taken from ED/brain subset)
    fit <- glmFit(DGE, dispersion = 0.3) 
    treat <- glmTreat(fit, coef = 2, lfc = 1) 
    summary(decideTests(treat, p = 0.05, adjust = "none"))
    
    #  
    results <- as.data.frame(topTags(treat, n = dim(DGE)[1]))
    sig_results <- results[results$PValue < 0.05,  ]
    
    write.csv(sig_results, file = "EdgeR_DEG_analysis/EdgeR_output_files/2.11.2021/26ed_v_8_I.brain.csv", row.names = FALSE)
    
    
## Output data ----
    

  