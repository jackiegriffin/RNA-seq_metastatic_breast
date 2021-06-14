# libraries ----
  library(clusterProfiler)
  library(magrittr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)

# adrenal data ----
  d_adrenal <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/adrenal.2.11.2021.csv")


# Format  ----
  geneList <- d_adrenal[,2]
  names(geneList) <- as.character(d_adrenal[,1])# feature 2: named vector
  geneList <- sort(geneList, decreasing = TRUE)# feature 3: decreasing order
  gene <- names(geneList)[abs(geneList) > 3] # set genelist to include those with FC > 2

# KEGG ----

  enrichkegg_adrenal <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05)

  enrichkegg_adrenal_plot <- dotplot(enrichkegg_adrenal, showCategory=5) 
    # ggtitle("KEGG Pathway Enrichment", subtitle = "Adrenal Metastases")
  plot(enrichkegg_adrenal_plot)
  
  
  enrichkegg_adrenal <- setReadable(enrichkegg_adrenal, 'org.Hs.eg.db', 'ENTREZID') 
  enrichkegg_adrenal_cnet_plot <- cnetplot(enrichkegg_adrenal, foldChange=geneList)  
    # ggtitle("KEGG Pathway Enrichment", subtitle = "Adrenal Metastases")
  plot(enrichkegg_adrenal_cnet_plot)


grid <- plot_grid(enrichKEGG_barplot, enrichkegg_adrenal_cnet_plot, ncol = 1) + ggtitle("Asdfs")


enrichKEGG_barplot <- barplot(enrichkegg_adrenal, showCategory=5)


plot(grid)
# cellular component enrichment (UPSET FOLDCHANGE PLOT) ----
  ego3_adrenal <- gseGO(geneList     = geneList,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "CC",
                      nPerm        = 1000,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      # pvalueCutoff = 0.05,
                      verbose      = FALSE)
  head(ego3_adrenal)
  upsetplot(ego3_adrenal)
  
  


