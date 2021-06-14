# libraries ----

  library(clusterProfiler)
  library(magrittr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(pathview)
  library(ggplot2)
  library(tibble)
  library(forcats)
  library(ggstance)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(rlang)
  library(tidyverse)

    # BiocManager::install('DOSE')
    # BiocManager::install('clusterProfiler')
    # BiocManager::install('tibble')
    # devtools::install_github("tidyverse/dplyr")

# adrenal data ----
  d_adrenal <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/Input_data/adrenal.2.11.2021.csv")


# Format  ----
  geneList <- d_adrenal[,2]
  names(geneList) <- as.character(d_adrenal[,1])# feature 2: named vector
  geneList <- sort(geneList, decreasing = TRUE)# feature 3: decreasing order
  gene <- names(geneList)[abs(geneList) > 3] # set genelist to include those with FC > 2

  
# enrichKEGG ----
  enrichkegg_adrenal <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05)
  enrichkegg_adrenal_readable <- setReadable(enrichkegg_adrenal, 'org.Hs.eg.db', 'ENTREZID') 
  head(enrichkegg_adrenal_readable)

   # splitplot 
  
      ggplot(enrichkegg_adrenal_readable, showCategory = 5, 
           aes(GeneRatio, fct_reorder(Description, GeneRatio))) + 
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=pvalue, size = Count)) +
        scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE), direction = -1) +
        scale_size_continuous(range=c(2, 10)) +
      theme_minimal() + 
      xlab("Gene Ratio (data/signature set)") +
      ylab(NULL) + 
      labs(caption = " ") +
        theme(axis.text.x = element_text(size = 18, colour = 'black')) +
        theme(axis.text.y = element_text(size = 22, color = 'black'))+
      theme(axis.title = element_text(size = 18)) +
      theme(title = element_text(size = 16)) +
      theme(plot.subtitle = element_text(size = 12)) +
      ggtitle("Pathway Enrichment in Adrenal Metastases", subtitle = "DEG Foldchange > 3") +
      theme(plot.caption = element_text(size = 10))
      
      
      
  # heatplot
      p2 <- heatplot(enrichkegg_adrenal_readable, foldChange=geneList)
      plot(p2) 
      
   # dotplot
      enrichkegg_adrenal_dotplot <- dotplot(enrichkegg_adrenal, showCategory=5) 
        # ggtitle("KEGG Pathway Enrichment", subtitle = "Adrenal Metastases")
      plot(enrichkegg_adrenal_dotplot)
    
  # cnetplot
      enrichkegg_adrenal <- setReadable(enrichkegg_adrenal, 'org.Hs.eg.db', 'ENTREZID') 
      enrichkegg_adrenal_cnet_plot <- cnetplot(enrichkegg_adrenal, foldChange=geneList)  
        # ggtitle("KEGG Pathway Enrichment", subtitle = "Adrenal Metastases")
      plot(enrichkegg_adrenal_cnet_plot)
  
  # barplot
      enrichKEGG_adrenal_barplot <- barplot(enrichkegg_adrenal, showCategory=5)
      plot(enrichKEGG_adrenal_barplot) + ylab("Gene Count") 
      
  # grid_plot
      grid <- plot_grid(enrichkegg_adrenal_dotplot, enrichkegg_adrenal_cnet_plot, ncol = 1) + ggtitle("Asdfs")
      plot(grid)


# gseGO CC UPSET plot FC ----
  ego3_adrenal <- gseGO(geneList = geneList,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "CC",
                    nPerm        = 1000,
                    minGSSize    = 100,
                    maxGSSize    = 500,
                    # pvalueCutoff = 0.05,
                    verbose      = FALSE)
  upsetplot(ego3_adrenal) + 
    ggtitle("Functional Enrichment for Cellular Component (gseGO)", subtitle = "Adrenal Metastases")
  
# pathview ----
  
  # hsa04151 = pi3Ki pathway
      hsa04151 <- pathview(gene.data  = geneList,
                           pathway.id = "hsa04151",
                           species    = "hsa",
                           limit      = list(gene=max(abs(geneList)), cpd=1),
                           out.suffix = "adrenal"
                           # kegg.native = FALSE)
      )
  
  # hsa04510 = focal adhesion 
      hsa04510 <- pathview(gene.data  = geneList,
                               pathway.id = "hsa04510",
                               species    = "hsa",
                               limit      = list(gene=max(abs(geneList)), cpd=1),
                               out.suffix = "adrenal",
                               kegg.native = FALSE)




      # hsa04512 =ECM recpetorinteractions
      
      hsa04512 <- pathview(gene.data  = geneList,
                           pathway.id = "hsa04512",
                           species    = "hsa",
                           limit      = list(gene=max(abs(geneList)), cpd=1),
                           out.suffix = "adrenal_ECM",
                           kegg.native = FALSE)





