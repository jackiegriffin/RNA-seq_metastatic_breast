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

# brain data ----
  d_brain <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/brain.2.11.2021.csv")

  

# Format  ----
  geneList <- d_brain[,2]
  names(geneList) <- as.character(d_brain[,1])
  geneList <- sort(geneList, decreasing = TRUE)
  gene <- names(geneList)[abs(geneList) > 3] 

  
# enrichKEGG ----
  enrichkegg_brain <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05)
  enrichkegg_brain_readable <- setReadable(enrichkegg_brain, 'org.Hs.eg.db', 'ENTREZID') 
  head(enrichkegg_brain_readable)
  
  # splitplot 
  
  ggplot(enrichkegg_brain_readable, showCategory = 8, 
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
    ggtitle("Pathway Enrichment in Brain Metastases", subtitle = "DEG Foldchange > 3") +
    theme(plot.caption = element_text(size = 10))
  
  
    # dotplot
    enrichkegg_brain_dotplot <- dotplot(enrichkegg_brain, showCategory=5) 
    # ggtitle("KEGG Pathway Enrichment", subtitle = "brain Metastases")
    plot(enrichkegg_brain_dotplot)
    
    # cnetplot
    enrichkegg_brain <- setReadable(enrichkegg_brain, 'org.Hs.eg.db', 'ENTREZID') 
    enrichkegg_brain_cnet_plot <- cnetplot(enrichkegg_brain_readable, foldChange=geneList)  
    # ggtitle("KEGG Pathway Enrichment", subtitle = "brain Metastases")
    plot(enrichkegg_brain_cnet_plot)
    
    # barplot
    enrichKEGG_brain_barplot <- barplot(enrichkegg_brain, showCategory=5)
    plot(enrichKEGG_brain_barplot) + ylab("Gene Count") 
    
    # grid_plot
    grid <- plot_grid(enrichkegg_brain_dotplot, enrichkegg_brain_cnet_plot, ncol = 1) + ggtitle("Asdfs")
    plot(grid)
    
    
# gseGO CC UPSET plot FC ----
    ego3_brain <- gseGO(geneList = geneList,
                          OrgDb        = org.Hs.eg.db,
                          ont          = "CC",
                          nPerm        = 1000,
                          minGSSize    = 100,
                          maxGSSize    = 500,
                          # pvalueCutoff = 0.05,
                          verbose      = FALSE)
    upsetplot(ego3_brain) + 
      ggtitle("Functional Enrichment for Cellular Component (gseGO)", subtitle = "Brain Metastases")
    
    
# pathview ----

    # hregulation of actin cytoskeleton = hsa04810  
      hsa04810 <- pathview(gene.data  = geneList,
                           pathway.id = "hsa04810",
                           species    = "hsa",
                           limit      = list(gene=max(abs(geneList)), cpd=1),
                           out.suffix = "brain_",
                           kegg.native = FALSE)
    
    # hsa04510 = focal adhesion 
      hsa04510 <- pathview(gene.data  = geneList,
                               pathway.id = "hsa04510",
                               species    = "hsa",
                               limit      = list(gene=max(abs(geneList)), cpd=1),
                               out.suffix = "brain",
                               kegg.native = FALSE)
          
        
      # Glioma = hsa05214 
        hsa05214 <- pathview(gene.data  = geneList,
                           pathway.id = "hsa05214",
                           species    = "hsa",
                           limit      = list(gene=max(abs(geneList)), cpd=1),
                           out.suffix = "brain")
                           # kegg.native = FALSE)
    


