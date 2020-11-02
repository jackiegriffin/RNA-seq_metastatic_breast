library(magrittr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

## ***************************** Pathway    ***********************************

## Upload data ----
d <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/gene.list_sig_complete_DEGs.csv")

## Format data ----
geneList <- d[,2]
names(geneList) <- as.character(d[,1])## feature 2: named vector
geneList <- sort(geneList, decreasing = TRUE)## feature 3: decreasing order
gene <- names(geneList)[abs(geneList) > 2] # set genelist to include those with FC > 2



## Enrichment plots ----

# DisGeNET(Janet et al. 2015) is an integrative and comprehensive resources 
# of gene-disease associations from several public data sources and the literature. 
# It contains gene-disease associations and snp-gene-disease associations.

# Dotplots 
genes_for_enrichDGN <- names(geneList)[abs(geneList) > 3] # enrichDGN
enrich_DGN <- enrichDGN(genes_for_enrichDGN)
p1 <- dotplot(enrich_DGN, showCategory=20) + ggtitle("enrichDGN, FC>3")

genes_for_enrichKEGG <- names(geneList)[abs(geneList) > 3] # enrichKEGG
enrich_KEGG <- enrichKEGG(gene = genes_for_enrichKEGG, organism='hsa',pvalueCutoff = 0.05)
p2 <- dotplot(enrich_KEGG, showCategory=20) + ggtitle("enrichKEGG, FC>3")
plot_grid(p1, p2, ncol=2)

# Gene network map 
enrich_DGN_cnet <- setReadable(enrich_DGN, 'org.Hs.eg.db', 'ENTREZID') # enrichDGN
p3 <- cnetplot(enrich_DGN_cnet, foldChange=geneList)

enrich_KEGG_cnet <- setReadable(enrich_KEGG, 'org.Hs.eg.db', 'ENTREZID') # enrichKEGG
p4 <- cnetplot(enrich_KEGG_cnet, foldChange=geneList)
plot_grid(p3, p4, ncol=2)



#MSigDb analysis ----
# 
# msigdbr_species()
# # retrieve all human data sets
# 
# m_df <- msigdbr(species = "Homo sapiens")
# head(m_df, 2) %>% as.data.frame
# 
# # select oncogenic C5 set as an example
# m_t2g <- msigdbr(species = "Homo sapiens", category = "C4") %>% 
#   dplyr::select(gs_name, entrez_gene)
# head(m_t2g)
# 
# 
# ## GSEA analysis ----
# em <- enricher(gene, TERM2GENE=m_t2g)
# em2 <- GSEA(geneList, TERM2GENE = m_t2g)
# em <- setReadable(em, org.Hs.eg.db, keyType = "ENTREZID") # HUGO gene names
# em2 <- setReadable(em2, org.Hs.eg.db, keyType = "ENTREZID") # HUGO gene names
# # head(em)
# # head(em2)
# # write.csv(em2, file = "em2_complete_sig.csv")
# # write.csv(em, file = "em_complete_sig.csv")
# 
# p3 <- dotplot(em, showCategory=30) + ggtitle("enriche_C4 ")
# p4 <- dotplot(em2, showCategory=30) + ggtitle("GSEA_C4")
# plot_grid(p3, p4, ncol=2)
# 
# 
# edox <- setReadable(em2, 'org.Hs.eg.db', 'ENTREZID')
# p2 <- cnetplot(edox, foldChange=geneList)
# plot(p2)


## Done ----

## ***************************** Pathway    ***********************************

##  Cluster profile by individual met ----
met_list <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/list_csv_trial1.csv", stringsAsFactors = FALSE, header = FALSE, row.names = 1)
met_list_matrix <- as.matrix(as.data.frame(met_list))
new_df<-strsplit(met_list_matrix,", ")
names(new_df) <- c("kidney.1", "Bone.4", "Bone.1", "Bone.3", "Adrenal.3", 
                   "Bone.2", "Adrenal.1", "Brain.4", "Bone.5", "Brain.1", 
                   "Brain.3", "Bone.6", "Brain.2", "Adrenal.2")



xx<- compareCluster(new_df, fun = "enrichDGN", organism="hsa", pvalueCutoff=0.05)

p5 <- emapplot(xx)
plot(p5)


