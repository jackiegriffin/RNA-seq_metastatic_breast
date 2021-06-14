# Load libraries ----
      library(biomaRt)
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
      library(viridis)
      library(ggsci)
      library(rlang)
      library(biomaRt)

      mart <- useMart('ENSEMBL_MART_ENSEMBL') #start
      mart <- useDataset('hsapiens_gene_ensembl', mart)


#######################  OBJECTIVE  ##########################################
#                                                                            #
#    PLATFORM : RNA-seq                                                      #
#    PIPELINE : FASTQC-STAR-Xenofilter-featurecounts                         #
#    INPUT (pipeline output) : feature_counts.txt                            #
#    Differential expression analysis (DEGs) using edgeR:                    #
#    OUTPUT* : Transcripts upregulated in organ, specific metastatic tumors  #
#         * See NOTES TO SELF at end of script                               #
#                                                                            #
##############################################################################

# Load data ----
    featureCounts <- read.table(file = "3. edgeR_upset/EdgeR_ud/Data_input/feature_counts.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
    metadata <- read.csv(file = "3. edgeR_upset/EdgeR_ud/Data_input/metadata_complete.csv", stringsAsFactors = FALSE, header = FALSE)
    
    # format df
      colnames(featureCounts)
      featureCounts_trim <- featureCounts[ -c(1:5) ]
      colnames(featureCounts_trim)
      colnames(featureCounts_trim)[1:21] <-c("primary_m1e2_4P", "adrenal_m1e2_22", "bone_m1e2_2_II", 
                                                 "primary_m2e2_5P", "adrenal_m2e2_12_I", "brain_m2e2_E", "bone_m2e2_5_I",
                                                 "primary_m3e2_1P", "adrenal_m3e2_2", "kidney_m3e2_12", "tibia_m3e2_C", "femur_m3e2_15_I",
                                                 "primary_m4ed_2P", "bone_m4ed_16",
                                                 "primary_m5ed_3P", "brain_m5ed_11_II",
                                                 "primary_m6ed_25", "bone_m6ed_10", "brain_m6ed_6",
                                                 "primary_m7ed_26", "brain_m7ed_8_I")

# Extrapolate BCV from ED brain data (n=3)  ----
          # Subset ED data
            ED_counts <- featureCounts_trim[c(15:17, 19:21)]
            colnames(ED_counts)

          # Asign groups
            groups_ED <- c(1,2,1,2,1,2) # 1 = primary 2 = met

          # DGEList
            DGE_ED_subset <- DGEList(counts = ED_counts, group = groups_ED, genes = row.names(ED_counts))
            dim(DGE_ED_subset)

          # Filter
            DGE_ED_subset <- DGE_ED_subset[rowSums(cpm(DGE_ED_subset) > 0.1) >= 6, keep.lib.sizes = FALSE]
            dim(DGE_ED_subset)

          # Normalize (TMM)
            DGE_ED_subset <- calcNormFactors(DGE_ED_subset, method = "TMM")
            DGE_ED_subset$samples

          # BCV
            DGE_ED_subset <- estimateCommonDisp(DGE_ED_subset, verbose = TRUE) # BCV = 0.48
    
# Apply BCV to matched primary and metastatic tumor case studies ----
      # Complete analysis (PCA and library visualization) ----
          # Assign groups
            colnames(featureCounts_trim)
            groups_complete <- c(1,2,2, # 1 = primary 2 = met
                                 1,2,2,2,
                                 1,2,2,2,2,
                                 1,2,
                                 1,2,
                                 1,2,2,
                                 1,2) 
            
            # DGEList
              DGE_complete <- DGEList(counts = featureCounts_trim, group = Groups_complete, genes = row.names(featureCounts_trim))
              dim(DGE_complete)
            
            # Filter
              DGE_complete_filt <- DGE_complete[rowSums(cpm(DGE_complete) > 0.1) >= 2, keep.lib.sizes = FALSE] 
              dim(DGE_complete_filt)
            
            # Normalize (TMM)
              DGE_complete_filt <- calcNormFactors(DGE_complete_filt, method = "TMM")
              DGE_complete_filt$samples
              
            # Plot counts & MDS plots
              cols <- as.numeric(DGE_complete_filt$samples$group)+2
              par(mar=c(4,4,4,3))
              barplot(colSums(DGE_complete_filt$counts), las=2, main="Counts per index",
                      col=cols, cex.names=0.5, cex.axis=0.8)
                      legend("topright", legend=c("Primary", "Metastases"), col=c(3,4), pch=15)
              cols2 <- DGE_complete_filt$samples$Groups_complete
              plotMDS(DGE_complete_filt, col=cols, main="MDS Plot")
                      legend("bottomright", legend=c("Primary", "metastases"), col=c(3,4), pch=15)
              
      # Individual case study analysis ----   
                      
          colnames(featureCounts_trim)
         
         # Subset df's            
            counts__4P_vs_22    <- featureCounts_trim[c(1,2)]    # E2 primary (4P) vs adrenal (22)
            counts__4P_vs_2.II  <- featureCounts_trim[c(1,3)]    # E2 primary (4P) vs bone (2_II)
            counts__5P_vs_12.I  <- featureCounts_trim[c(4,5)]    # E2 primary (5P) vs adrenal (12_I)
            counts__5P_vs_E     <- featureCounts_trim[c(4,6)]    # E2 primary (5P) vs brain (E)
            counts__5P_vs_5.I   <- featureCounts_trim[c(4,7)]    # E2 primary (5P) vs bone (5.I)
            counts__1P_vs_2     <- featureCounts_trim[c(8,9)]    # E2 primary (1P) vs adrenal (2)
            counts__1P_vs_12    <- featureCounts_trim[c(8,10)]   # E2 primary (1P) vs kidney (12)
            counts__1P_vs_C     <- featureCounts_trim[c(8,11)]   # E2 primary (1P) vs bone (C)
            counts__1P_vs_15.I  <- featureCounts_trim[c(8,12)]   # E2 primary (1P) vs bone (15.I)
            counts__2P_vs_16    <- featureCounts_trim[c(13,14)]  # ED primary (2P) vs bone (16)
            counts__3P_vs_11.II <- featureCounts_trim[c(15,16)]  # ED primary (3P) vs bone (11.II)
            counts__25_vs_10    <- featureCounts_trim[c(17,18)]  # ED primary (25) vs bone (10)
            counts__25_vs_6     <- featureCounts_trim[c(17,19)]  # ED primary (25) vs brain (6)
            counts__26_vs_8.I   <- featureCounts_trim[c(20,21)]  # ED primary (26) vs brain (8.I)
            
          # Assign groups 
            groups_case_studies <- c(1,2) # 1 = primary 2 = met
  
          # DGE models
            # counts__4P_vs_22 ----
  
                        # DGEList
                          counts__4P_vs_22<- DGEList(counts = counts__4P_vs_22, group = groups_case_studies, genes = row.names(counts))
                          dim(DGE__4P_vs_22)
                        
                        # Filter 
                          counts__4P_vs_22 <- counts__4P_vs_22[rowSums(cpm(counts__4P_vs_22) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__4P_vs_22)
                        
                        # Normalize (TMM)
                          counts__4P_vs_22 <- calcNormFactors(counts__4P_vs_22, method = "TMM")
                          counts__4P_vs_22$samples
                          plotMD(cpm(counts__4P_vs_22, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__4P_vs_22, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                        
                        # Results
                          results__4P_vs_22 <- as.data.frame(topTags(treat, n = dim(counts__4P_vs_22)[1]))
                          significant__results__4P_vs_22 <- results__4P_vs_22[results__4P_vs_22$PValue < 0.05,  ]                   # Subset significant p<0.05
                          significant__results__4P_vs_22 <- tibble::rownames_to_column(significant__results__4P_vs_22, "row_names") # Apply row.names_to_column
                          colnames(significant__results__4P_vs_22)[1] <-c("Geneid")                                                 # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc_4P_vs_22 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results__4P_vs_22$Geneid) 
                            significant__results__4P_vs_22_hgnc <- merge(significant__results__4P_vs_22, annotate_hgnc_4P_vs_22, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results__4P_vs_22_hgnc$sample <- rep("Adrenal.22", nrow(significant__results__4P_vs_22_hgnc)) # Add sample id info into column
                            head(significant__results__4P_vs_22_hgnc) # check 
                            
                            # Write .csv file 
                              write.csv(significant__results__4P_vs_22_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results__4P_vs_22_hgnc.csv", row.names = TRUE)
                          
                          
            # counts__4P_vs_2.II ----
                          
                        # DGEList
                          counts__4P_vs_2.II<- DGEList(counts = counts__4P_vs_2.II, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__4P_vs_2.II)
                          
                        # Filter 
                          counts__4P_vs_2.II <- counts__4P_vs_2.II[rowSums(cpm(counts__4P_vs_2.II) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__4P_vs_2.II)
                          
                        # Normalize (TMM)
                          counts__4P_vs_2.II <- calcNormFactors(counts__4P_vs_2.II, method = "TMM")
                          counts__4P_vs_2.II$samples
                          plotMD(cpm(counts__4P_vs_2.II, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__4P_vs_2.II, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results__4P_vs_2.II <- as.data.frame(topTags(treat, n = dim(counts__4P_vs_2.II)[1]))
                          significant__results__4P_vs_2.II <- results__4P_vs_2.II[results__4P_vs_2.II$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results__4P_vs_2.II <- tibble::rownames_to_column(significant__results__4P_vs_2.II, "row_names") # Apply row.names_to_column
                          colnames(significant__results__4P_vs_2.II)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__4P_vs_2.II <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results__4P_vs_2.II$Geneid) 
                            significant__results__4P_vs_2.II_hgnc <- merge(significant__results__4P_vs_2.II, annotate_hgnc__4P_vs_2.II, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results__4P_vs_2.II_hgnc$sample <- rep("Bone.2_II", nrow(significant__results__4P_vs_2.II_hgnc)) # Add sample id info into column
                            head(significant__results__4P_vs_2.II_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results__4P_vs_2.II_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results__4P_vs_2.II_hgnc.csv", row.names = FALSE)
                          
            # counts__5P_vs_12.I ----
                          
                        # DGEList
                          counts__5P_vs_12.I<- DGEList(counts = counts__5P_vs_12.I, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__5P_vs_12.I)
                          
                        # Filter 
                          counts__5P_vs_12.I <- counts__5P_vs_12.I[rowSums(cpm(counts__5P_vs_12.I) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__5P_vs_12.I)
                          
                        # Normalize (TMM)
                          counts__5P_vs_12.I <- calcNormFactors(counts__5P_vs_12.I, method = "TMM")
                          counts__5P_vs_12.I$samples
                          plotMD(cpm(counts__5P_vs_12.I, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__5P_vs_12.I, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results 
                          results__5P_vs_12.I <- as.data.frame(topTags(treat, n = dim(counts__5P_vs_12.I)[1]))
                          significant__results__5P_vs_12.I <- results__5P_vs_12.I[results__5P_vs_12.I$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results__5P_vs_12.I <- tibble::rownames_to_column(significant__results__5P_vs_12.I, "row_names") # Apply row.names_to_column
                          colnames(significant__results__5P_vs_12.I)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__5P_vs_12.I <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results__5P_vs_12.I$Geneid) 
                            significant__results__5P_vs_12.I_hgnc <- merge(significant__results__5P_vs_12.I, annotate_hgnc__5P_vs_12.I, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results__5P_vs_12.I_hgnc$sample <- rep("Adrenal.12_II", nrow(significant__results__5P_vs_12.I_hgnc)) # Add sample id info into column
                            head(significant__results__5P_vs_12.I_hgnc) # check
                          
                            # Write .csv file 
                              write.csv(significant__results__5P_vs_12.I_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results__5P_vs_12.I_hgnc.csv", row.names = FALSE)
                          
                          
            # counts__5P_vs_E ----
                          
                        # DGEList
                          counts__5P_vs_E<- DGEList(counts = counts__5P_vs_E, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__5P_vs_E)
                          
                        # Filter 
                          counts__5P_vs_E <- counts__5P_vs_E[rowSums(cpm(counts__5P_vs_E) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__5P_vs_E)
                          
                        # Normalize (TMM)
                          counts__5P_vs_E <- calcNormFactors(counts__5P_vs_E, method = "TMM")
                          counts__5P_vs_E$samples
                          plotMD(cpm(counts__5P_vs_E, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__5P_vs_E, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____5P_vs_E <- as.data.frame(topTags(treat, n = dim(counts__5P_vs_E)[1]))
                          significant__results__5P_vs_E <- results____5P_vs_E[results____5P_vs_E$PValue < 0.05,  ]                # Subset significant p<0.05
                          significant__results__5P_vs_E <- tibble::rownames_to_column(significant__results__5P_vs_E, "row_names") # Apply row.names_to_column
                          colnames(significant__results__5P_vs_E)[1] <-c("Geneid")                                                # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__5P_vs_E <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results__5P_vs_E$Geneid) 
                            significant__results__5P_vs_E_hgnc <- merge(significant__results__5P_vs_E, annotate_hgnc__5P_vs_E, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results__5P_vs_E_hgnc$sample <- rep("Brain.E", nrow(significant__results__5P_vs_E_hgnc)) # Add sample id info into column
                            head(significant__results__5P_vs_E_hgnc) # check
                            
                            # Write .csv file 
                              write.csv(significant__results__5P_vs_E_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results__5P_vs_E_hgnc.csv", row.names = FALSE)
                          
                          
            # counts__5P_vs_5.I ----
                          
                        # DGEList
                          counts__5P_vs_5.I<- DGEList(counts = counts__5P_vs_5.I, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__5P_vs_5.I)
                          
                        # Filter 
                          counts__5P_vs_5.I <- counts__5P_vs_5.I[rowSums(cpm(counts__5P_vs_5.I) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__5P_vs_5.I)
                          
                        # Normalize (TMM)
                          counts__5P_vs_5.I <- calcNormFactors(counts__5P_vs_5.I, method = "TMM")
                          counts__5P_vs_5.I$samples
                          plotMD(cpm(counts__5P_vs_5.I, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__5P_vs_5.I, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results 
                          results__5P_vs_5.I <- as.data.frame(topTags(treat, n = dim(counts__5P_vs_5.I)[1]))
                          significant__results__5P_vs_5.I <- results__5P_vs_5.I[results__5P_vs_5.I$PValue < 0.05,  ]                  # Subset significant p<0.05
                          significant__results__5P_vs_5.I <- tibble::rownames_to_column(significant__results__5P_vs_5.I, "row_names") # Apply row.names_to_column
                          colnames(significant__results__5P_vs_5.I)[1] <-c("Geneid")                                                  # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__5P_vs_5.I <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results__5P_vs_5.I$Geneid) 
                            significant__results__5P_vs_5.I_hgnc <- merge(significant__results__5P_vs_5.I, annotate_hgnc__5P_vs_5.I, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results__5P_vs_5.I_hgnc$sample <- rep("Bone.5_I", nrow(significant__results__5P_vs_5.I_hgnc)) # Add sample id info into column
                            head(significant__results__5P_vs_5.I_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results__5P_vs_5.I_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results__5P_vs_5.I_hgnc.csv", row.names = FALSE)
                          
            # counts__1P_vs_2 ----
                          
                        # DGEList
                          counts__1P_vs_2<- DGEList(counts = counts__1P_vs_2, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__1P_vs_2)
                          
                        # Filter 
                          counts__1P_vs_2 <- counts__1P_vs_2[rowSums(cpm(counts__1P_vs_2) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__1P_vs_2)
                          
                        # Normalize (TMM)
                          counts__1P_vs_2 <- calcNormFactors(counts__1P_vs_2, method = "TMM")
                          counts__1P_vs_2$samples
                          plotMD(cpm(counts__1P_vs_2, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__1P_vs_2, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____1P_vs_2 <- as.data.frame(topTags(treat, n = dim(counts__1P_vs_2)[1]))
                          significant__results____1P_vs_2 <- results____1P_vs_2[results____1P_vs_2$PValue < 0.05,  ]                  # Subset significant p<0.05
                          significant__results____1P_vs_2 <- tibble::rownames_to_column(significant__results____1P_vs_2, "row_names") # Apply row.names_to_column
                          colnames(significant__results____1P_vs_2)[1] <-c("Geneid")                                                  # Re-name row.names column to Geneid
                          head(significant__results____1P_vs_2)
                         
                           # Add column with hugo gene names to df
                            annotate_hgnc__1P_vs_2 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____1P_vs_2$Geneid) 
                            significant__results____1P_vs_2_hgnc <- merge(significant__results____1P_vs_2, annotate_hgnc__1P_vs_2, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____1P_vs_2_hgnc$sample <- rep("Adrenal.2", nrow(significant__results____1P_vs_2_hgnc)) # Add sample id info into column
                            head(significant__results____1P_vs_2_hgnc) # check

                            # Write .csv file 
                              write.csv(significant__results____1P_vs_2_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____1P_vs_2_hgnc.csv", row.names = FALSE)
                          
                          
            # counts__1P_vs_12 ----
                          
                        # DGEList
                          counts__1P_vs_12<- DGEList(counts = counts__1P_vs_12, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__1P_vs_12)
                          
                        # Filter 
                          counts__1P_vs_12 <- counts__1P_vs_12[rowSums(cpm(counts__1P_vs_12) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__1P_vs_12)
                          
                        # Normalize (TMM)
                          counts__1P_vs_12 <- calcNormFactors(counts__1P_vs_12, method = "TMM")
                          counts__1P_vs_12$samples
                          plotMD(cpm(counts__1P_vs_12, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__1P_vs_12, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____1P_vs_12 <- as.data.frame(topTags(treat, n = dim(counts__1P_vs_12)[1]))
                          significant__results____1P_vs_12 <- results____1P_vs_12[results____1P_vs_12$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results____1P_vs_12 <- tibble::rownames_to_column(significant__results____1P_vs_12, "row_names") # Apply row.names_to_column
                          colnames(significant__results____1P_vs_12)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__1P_vs_12 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____1P_vs_12$Geneid) 
                            significant__results____1P_vs_12_hgnc <- merge(significant__results____1P_vs_2, annotate_hgnc__1P_vs_12, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____1P_vs_12_hgnc$sample <- rep("Kidney.12", nrow(significant__results____1P_vs_12_hgnc)) # Add sample id info into column
                            head(significant__results____1P_vs_12_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____1P_vs_12_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____1P_vs_12_hgnc.csv", row.names = FALSE)
                          
            # counts__1P_vs_C ----
                          
                        # DGEList
                          counts__1P_vs_C<- DGEList(counts = counts__1P_vs_C, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__1P_vs_C)
                          
                        # Filter 
                          counts__1P_vs_C <- counts__1P_vs_C[rowSums(cpm(counts__1P_vs_C) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__1P_vs_C)
                          
                        # Normalize (TMM)
                          counts__1P_vs_C <- calcNormFactors(counts__1P_vs_C, method = "TMM")
                          counts__1P_vs_C$samples
                          plotMD(cpm(counts__1P_vs_C, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__1P_vs_C, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____1P_vs_C <- as.data.frame(topTags(treat, n = dim(counts__1P_vs_C)[1]))
                          significant__results____1P_vs_C <- results____1P_vs_C[results____1P_vs_C$PValue < 0.05,  ]                   # Subset significant p<0.05
                          significant__results____1P_vs_C <- tibble::rownames_to_column(significant__results____1P_vs_C, "row_names")  # Apply row.names_to_column
                          colnames(significant__results____1P_vs_C)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__1P_vs_C <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____1P_vs_C$Geneid) 
                            significant__results____1P_vs_C_hgnc <- merge(significant__results____1P_vs_C, annotate_hgnc__1P_vs_C, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____1P_vs_C_hgnc$sample <- rep("Bone.C", nrow(significant__results____1P_vs_C_hgnc)) # Add sample id info into column
                            head(significant__results____1P_vs_C_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____1P_vs_C_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____1P_vs_C_hgnc.csv", row.names = FALSE)
                              
            # counts__1P_vs_15.I ----
                          
                        # DGEList
                          counts__1P_vs_15.I<- DGEList(counts = counts__1P_vs_15.I, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__1P_vs_15.I)
                          
                        # Filter 
                          counts__1P_vs_15.I <- counts__1P_vs_15.I[rowSums(cpm(counts__1P_vs_15.I) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__1P_vs_15.I)
                          
                        # Normalize (TMM)
                          counts__1P_vs_15.I <- calcNormFactors(counts__1P_vs_15.I, method = "TMM")
                          counts__1P_vs_15.I$samples
                          plotMD(cpm(counts__1P_vs_15.I, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__1P_vs_15.I, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____1P_vs_15 <- as.data.frame(topTags(treat, n = dim(counts__1P_vs_15.I)[1]))
                          significant__results____1P_vs_15 <- results____1P_vs_15[results____1P_vs_15$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results____1P_vs_15 <- tibble::rownames_to_column(significant__results____1P_vs_15, "row_names") # Apply row.names_to_column
                          colnames(significant__results____1P_vs_15)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__1P_vs_15.I <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____1P_vs_15$Geneid) 
                            significant__results____1P_vs_15.I_hgnc <- merge(significant__results____1P_vs_15, annotate_hgnc__1P_vs_15.I, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____1P_vs_15.I_hgnc$sample <- rep("Bone.15.I", nrow(significant__results____1P_vs_15.I_hgnc)) # Add sample id info into column
                            head(significant__results____1P_vs_15.I_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____1P_vs_15.I_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____1P_vs_15.I_hgnc.csv", row.names = FALSE)
                          
                          
                          
            # counts__2P_vs_16 ----
                          
                        # DGEList
                          counts__2P_vs_16<- DGEList(counts = counts__2P_vs_16, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__2P_vs_16)
                          
                        # Filter 
                          counts__2P_vs_16 <- counts__2P_vs_16[rowSums(cpm(counts__2P_vs_16) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__2P_vs_16)
                          
                        # Normalize (TMM)
                          counts__2P_vs_16 <- calcNormFactors(counts__2P_vs_16, method = "TMM")
                          counts__2P_vs_16$samples
                          plotMD(cpm(counts__2P_vs_16, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__2P_vs_16, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____2P_vs_16 <- as.data.frame(topTags(treat, n = dim(counts__2P_vs_16)[1]))
                          significant__results____2P_vs_16 <- results____2P_vs_16[results____2P_vs_16$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results____2P_vs_16 <- tibble::rownames_to_column(significant__results____2P_vs_16, "row_names") # Apply row.names_to_column
                          colnames(significant__results____2P_vs_16)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__2P_vs_16 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____2P_vs_16$Geneid) 
                            significant__results____2P_vs_16_hgnc <- merge(significant__results____2P_vs_16, annotate_hgnc__2P_vs_16, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____2P_vs_16_hgnc$sample <- rep("Bone.16", nrow(significant__results____2P_vs_16_hgnc)) # Add sample id info into column
                            head(significant__results____2P_vs_16_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____2P_vs_16_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____2P_vs_16_hgnc.csv", row.names = FALSE)
                          
            # counts__3P_vs_11.II ----
                          
                        # DGEList
                          counts__3P_vs_11.II<- DGEList(counts = counts__3P_vs_11.II, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__3P_vs_11.II)
                          
                        # Filter 
                          counts__3P_vs_11.II <- counts__3P_vs_11.II[rowSums(cpm(counts__3P_vs_11.II) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__3P_vs_11.II)
                          
                        # Normalize (TMM)
                          counts__3P_vs_11.II <- calcNormFactors(counts__3P_vs_11.II, method = "TMM")
                          counts__3P_vs_11.II$samples
                          plotMD(cpm(counts__3P_vs_11.II, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__3P_vs_11.II, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____3P_vs_11.II <- as.data.frame(topTags(treat, n = dim(counts__3P_vs_11.II)[1]))
                          significant__results____3P_vs_11.II <- results____3P_vs_11.II[results____3P_vs_11.II$PValue < 0.05,  ]              # Subset significant p<0.05
                          significant__results____3P_vs_11.II <- tibble::rownames_to_column(significant__results____3P_vs_11.II, "row_names") # Apply row.names_to_column
                          colnames(significant__results____3P_vs_11.II)[1] <-c("Geneid")                                                      # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__3P_vs_11.II <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____3P_vs_11.II$Geneid) 
                            significant__results____3P_vs_11.II_hgnc <- merge(significant__results____3P_vs_11.II, annotate_hgnc__3P_vs_11.II, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____3P_vs_11.II_hgnc$sample <- rep("Bone.11.II", nrow(significant__results____3P_vs_11.II_hgnc)) # Add sample id info into column
                            head(significant__results____3P_vs_11.II_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____3P_vs_11.II_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____3P_vs_11.II_hgnc.csv", row.names = FALSE)
                          
            # counts__25_vs_10 ----
                          
                        # DGEList
                          counts__25_vs_10<- DGEList(counts = counts__25_vs_10, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__25_vs_10)
                          
                        # Filter 
                          counts__25_vs_10 <- counts__25_vs_10[rowSums(cpm(counts__25_vs_10) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__25_vs_10)
                          
                        # Normalize (TMM)
                          counts__25_vs_10 <- calcNormFactors(counts__25_vs_10, method = "TMM")
                          counts__25_vs_10$samples
                          plotMD(cpm(counts__25_vs_10, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__25_vs_10, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____25_vs_10 <- as.data.frame(topTags(treat, n = dim(counts__25_vs_10)[1]))
                          significant__results____25_vs_10 <- results____25_vs_10[results____25_vs_10$PValue < 0.05,  ]                 # Subset significant p<0.05
                          significant__results____25_vs_10 <- tibble::rownames_to_column(significant__results____25_vs_10, "row_names") # Apply row.names_to_column
                          colnames(significant__results____25_vs_10)[1] <-c("Geneid")                                                   # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__25_vs_10 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____25_vs_10$Geneid) 
                            significant__results____25_vs_10_hgnc <- merge(significant__results____25_vs_10, annotate_hgnc__25_vs_10, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____25_vs_10_hgnc$sample <- rep("Bone.10", nrow(significant__results____25_vs_10_hgnc)) # Add sample id info into column
                            head(significant__results____25_vs_10_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____25_vs_10_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____25_vs_10_hgnc.csv", row.names = FALSE)
            # counts__25_vs_6 ----
                          
                        # DGEList
                          counts__25_vs_6<- DGEList(counts = counts__25_vs_6, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__25_vs_6)
                          
                        # Filter 
                          counts__25_vs_6 <- counts__25_vs_6[rowSums(cpm(counts__25_vs_6) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__25_vs_6)
                          
                        # Normalize (TMM)
                          counts__25_vs_6 <- calcNormFactors(counts__25_vs_6, method = "TMM")
                          counts__25_vs_6$samples
                          plotMD(cpm(counts__25_vs_6, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__25_vs_6, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____25_vs_6 <- as.data.frame(topTags(treat, n = dim(counts__25_vs_6)[1]))
                          significant__results____25_vs_6 <- results____25_vs_6[results____25_vs_6$PValue < 0.05,  ]                  # Subset significant p<0.05
                          significant__results____25_vs_6 <- tibble::rownames_to_column(significant__results____25_vs_6, "row_names") # Apply row.names_to_column
                          colnames(significant__results____25_vs_6)[1] <-c("Geneid")                                                  # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__25_vs_6 <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____25_vs_6$Geneid) 
                            significant__results____25_vs_6_hgnc <- merge(significant__results____25_vs_6, annotate_hgnc__25_vs_6, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____25_vs_6_hgnc$sample <- rep("Brain.6", nrow(significant__results____25_vs_6_hgnc)) # Add sample id info into column
                            head(significant__results____25_vs_6_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____25_vs_6_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____25_vs_6_hgnc.csv", row.names = FALSE)
            # counts__26_vs_8.I ----
                          
                        # DGEList
                          counts__26_vs_8.I<- DGEList(counts = counts__26_vs_8.I, group = groups_case_studies, genes = row.names(counts))
                          dim(counts__26_vs_8.I)
                          
                        # Filter 
                          counts__26_vs_8.I <- counts__26_vs_8.I[rowSums(cpm(counts__26_vs_8.I) > .1) >= 2, keep.lib.sizes = FALSE] 
                          dim(counts__26_vs_8.I)
                          
                        # Normalize (TMM)
                          counts__26_vs_8.I <- calcNormFactors(counts__26_vs_8.I, method = "TMM")
                          counts__26_vs_8.I$samples
                          plotMD(cpm(counts__26_vs_8.I, log = TRUE), column = 1) # check performance
                          abline(h=0, col="red", lty=2, lwd=2)
                          
                        # glmFit (BCV taken from ED/brain subset)
                          fit <- glmFit(counts__26_vs_8.I, dispersion = 0.3) 
                          treat <- glmTreat(fit, coef = 2, lfc = 1) 
                          summary(decideTests(treat, p = 0.05, adjust = "none"))
                          
                        # Results
                          results____26_vs_8.I <- as.data.frame(topTags(treat, n = dim(counts__26_vs_8.I)[1]))
                          significant__results____26_vs_8.I <- results____26_vs_8.I[results____26_vs_8.I$PValue < 0.05,  ]                # Subset significant p<0.05
                          significant__results____26_vs_8.I <- tibble::rownames_to_column(significant__results____26_vs_8.I, "row_names") # Apply row.names_to_column
                          colnames(significant__results____26_vs_8.I)[1] <-c("Geneid")                                                    # Re-name row.names column to Geneid

                          # Add column with hugo gene names to df
                            annotate_hgnc__26_vs_8.I <- getBM(mart=mart, attributes = c('entrezgene_id', 'hgnc_symbol'), filter = 'entrezgene_id', values = significant__results____26_vs_8.I$Geneid) 
                            significant__results____26_vs_8.I_hgnc <- merge(significant__results____26_vs_8.I, annotate_hgnc__26_vs_8.I, by.x = "Geneid", by.y = "entrezgene_id")
                            significant__results____26_vs_8.I_hgnc$sample <- rep("Brain.8.I", nrow(significant__results____26_vs_8.I_hgnc)) # Add sample id info into column
                            head(significant__results____26_vs_8.I_hgnc) # check 
                          
                            # Write .csv file 
                              write.csv(significant__results____26_vs_8.I_hgnc, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/significant__results____26_vs_8.I_hgnc.csv", row.names = FALSE)
                          
                          
                          
                          

# Complete DEG List ----
  complete_DEGs <- rbind(significant__results__4P_vs_22_hgnc, significant__results__4P_vs_2.II_hgnc, significant__results__5P_vs_12.I_hgnc, significant__results__5P_vs_E_hgnc, significant__results__5P_vs_5.I_hgnc,
                         significant__results____1P_vs_2_hgnc, significant__results____1P_vs_12_hgnc, significant__results____1P_vs_C_hgnc, significant__results____1P_vs_15.I_hgnc, 
                         significant__results____2P_vs_16_hgnc, significant__results____3P_vs_11.II_hgnc, significant__results____25_vs_10_hgnc, significant__results____25_vs_6_hgnc, 
                         significant__results____26_vs_8.I_hgnc)
                              
                              write.csv(complete_DEGs, file = "3. edgeR_upset/EdgeR_ud/edgeR_DEG_output/complete_DEGs.csv", row.names = FALSE)
                              
                              
                            
                            


