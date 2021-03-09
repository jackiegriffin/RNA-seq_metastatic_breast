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


# bone data ----
  d_bone <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/bone.2.11.2021.csv")

# Format  ----
  geneList <- d_bone[,2]
  names(geneList) <- as.character(d_bone[,1])# feature 2: named vector
  geneList <- sort(geneList, decreasing = TRUE)# feature 3: decreasing order
  gene <- names(geneList)[abs(geneList) > 3] # set genelist to include those with FC > 2


# enrichKEGG ----
  enrichkegg_bone <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05)

  
  enrichkegg_bone_readable <- setReadable(enrichkegg_bone, 'org.Hs.eg.db', 'ENTREZID') 
  head(enrichkegg_bone_readable)
  
  # splitplot 
  
  ggplot(enrichkegg_bone_readable, showCategory = 8, 
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
    theme(axis.text.y = element_text(size = 28, color = 'black'))+
    theme(axis.title = element_text(size = 18)) +
    theme(title = element_text(size = 16)) +
    theme(plot.subtitle = element_text(size = 12)) +
    ggtitle("Pathway Enrichment in Bone Metastases", subtitle = "DEG Foldchange > 3") +
    theme(plot.caption = element_text(size = 10))
  
  
  
  
  # dotplot
    enrichkegg_bone_dotplot <- dotplot(enrichkegg_bone_readable, showCategory=5) 
  # ggtitle("KEGG Pathway Enrichment", subtitle = "bone Metastases")
    plot(enrichkegg_bone_dotplot)
  
  # cnetplot
    enrichkegg_bone <- setReadable(enrichkegg_bone, 'org.Hs.eg.db', 'ENTREZID') 
    enrichkegg_bone_cnet_plot <- cnetplot(enrichkegg_bone_readable, foldChange=geneList)  
  # ggtitle("KEGG Pathway Enrichment", subtitle = "bone Metastases")
    plot(enrichkegg_bone_cnet_plot)
  
  # barplot
    enrichKEGG_bone_barplot <- barplot(enrichkegg_bone, showCategory=5)
    plot(enrichKEGG_bone_barplot) + ylab("Gene Count") 
  
  # grid_plot
    grid <- plot_grid(enrichkegg_bone_dotplot, enrichkegg_bone_cnet_plot, ncol = 1) + ggtitle("Asdfs")
    plot(grid)


# gseGO CC UPSET plot FC ----
  ego3_bone <- gseGO(geneList = geneList,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "CC",
                      nPerm        = 1000,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      # pvalueCutoff = 0.05,
                      verbose      = FALSE)
  upsetplot(ego3_bone) + 
    ggtitle("Functional Enrichment for Cellular Component (gseGO)", subtitle = "bone Metastases")


# pathview ----

  # Pathways of neurodegeneration - multiple diseases - Homo sapiens (human) = hsa05022
    hsa05022 <- pathview(gene.data  = geneList,
                       pathway.id = "hsa05022",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneList)), cpd=1),
                       out.suffix = "bone_")

  # Amyotrophic lateral sclerosis - Homo sapiens (human) = 	hsa05014    
    hsa05014 <- pathview(gene.data  = geneList,
                       pathway.id = "hsa05014",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneList)), cpd=1),
                       out.suffix = "bone")





