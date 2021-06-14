# libraries ----
  library(magrittr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)

# DEGs ----
  gene_fc_complete <- read.csv("3. edgeR_upset/pathway/complete_DEG_FC.csv")

# Format  ----
  geneList <- gene_fc_complete[,2]
  names(geneList) <- as.character(gene_fc_complete[,1])# feature 2: named vector
  geneList <- sort(geneList, decreasing = TRUE)# feature 3: decreasing order
  gene <- names(geneList)[abs(geneList) > 3] # set genelist to include those with FC > 2

# GSEA, KEGG, disgenet ---- 
    gsea <- gseNCG(geneList, nPerm=10000)
    gsea_plot <- dotplot(gsea, showCategory=20) + ggtitle("Functional enrichment by cancer type (GSEA)")  
    plot(gsea_plot)
    enrich_KEGG <- enrichKEGG(gene = gene, organism='hsa',pvalueCutoff = 0.05)
    kegg_plot <- dotplot(enrich_KEGG, showCategory=8) + ggtitle("KEGG")
    
    enrich_DGN <- enrichDGN(gene)
    disgenet_plot <- dotplot(enrich_DGN, showCategory=8) + ggtitle("DisGeNET")
    
    plot_grid(gsea_plot, kegg_plot, disgenet_plot)
    ridgeplot(gsea) + ggtitle("GSEA directionality")
   
# pubmed ----
    terms <- enrich_KEGG$Description[1:8]
    pubmed_KEGG <- pmcplot(terms, 2010:2020) + ggtitle("Pubmed trends of enriched terms")
    plot(pubmed_KEGG)

    
# Gene network maps ----
    gsea_cnet <- setReadable(gsea, 'org.Hs.eg.db', 'ENTREZID') # enrichDGN
    gse_network <- cnetplot(gsea_cnet, foldChange=geneList)
    plot(gse_network)

    enrich_KEGG_cnet <- setReadable(enrich_KEGG, 'org.Hs.eg.db', 'ENTREZID') # enrichKEGG
    kegg_network <- cnetplot(enrich_KEGG_cnet, foldChange=geneList) + ggtitle("KEGG network")
    plot(kegg_network)


    enrich_disgenet_cnet <- setReadable(enrich_DGN, 'org.Hs.eg.db', 'ENTREZID') # enrichKEGG
    disgenet_network <- cnetplot(enrich_disgenet_cnet, foldChange=geneList) + ggtitle("Disgenet network")
    plot(disgenet_network)



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


# Format  ----
geneList <- gene_fc_complete[,2]
names(geneList) <- as.character(gene_fc_complete[,1])# feature 2: named vector
geneList <- sort(geneList, decreasing = TRUE)# feature 3: decreasing order
gene <- names(geneList)[abs(geneList) > 3] # set genelist to include those with FC > 2


# Network cancer gene enrichment ---- 
gsea <- gseNCG(geneList, nPerm=10000)

# ridgeplot gseNCG
NCG_ridgeplot <- ridgeplot(gsea, showCategory = 5) + 
  ggtitle("Functional Enrichment for Cancer Type (gseNCG)") + 
  xlab("Expression Fold Change")
plot(NCG_ridgeplot)


# gseGO CC----
ego3_complete <- gseGO(geneList     = geneList,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "CC",
                       nPerm        = 1000,
                       minGSSize    = 100,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE) 
# Upset plot   
upsetplot(ego3_complete) + 
  ggtitle("Functional Enrichment for Cellular Component (gseGO)", subtitle = "Complete DEG's")


# pathview ----

# hsa05224 = breast cancer

hsa05224 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05224",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1),
                     out.suffix = "complete3_DEGs_breast_cancer",
                     kegg.native = FALSE)




