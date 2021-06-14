# libraries ----
  library(magrittr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(enrichplot)

# DEGs ----
  gene_fc_complete <- read.csv("3. edgeR_upset/pathway/complete_DEG_FC.csv")

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
                
        
        
        
        
        
        