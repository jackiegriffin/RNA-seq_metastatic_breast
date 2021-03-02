# Libraries ----
library("pathview")

#   KEGG Pathway Maps ----
#       breast cancer = hsa05224
#       PI3K-Akt signaling pathway = hsa04151
#       Transcriptional misregulation in cancer = hsa05202
#       MicroRNAs in cancer = hsa05206
#      
##    BONE ONLY --------------------------------------------------------
  bone <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/bone.2.11.2021.csv")
  geneListbone <- bone[,2]
  names(geneListbone) <- as.character(bone[,1])
  geneListbone <- sort(geneListbone, decreasing = TRUE)
  gene <- names(geneListbone)[abs(geneListbone) > 1] 
  hsa05014 <- pathview(gene.data  = geneListbone,
                              pathway.id = "hsa05014",
                              species    = "hsa",
                              limit      = list(gene=max(abs(geneListbone)), cpd=1),
                              out.suffix = "bone")
  
  ##    BRAIN ONLY --------------------------------------------------------
  # brain <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/brain.2.11.2021.csv")
  # geneListbrain <- brain[,2]
  # names(geneListbrain) <- as.character(brain[,1])
  # geneListbrain <- sort(geneListbrain, decreasing = TRUE)
  # gene <- names(geneListbrain)[abs(geneListbrain) > 1] 
  # hsa05214 <- pathview(gene.data  = geneListbrain,
  #                           pathway.id = "hsa05214",
  #                           species    = "hsa",
  #                           limit      = list(gene=max(abs(geneListbrain)), cpd=1),
  #                           out.suffix = "brain",
  #                      kegg.native = FALSE)
  #                      # kegg.native = FALSE)
  
  # ##    ADRENAL ONLY --------------------------------------------------------
  adrenal <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/adrenal.2.11.2021.csv")
  geneListadrenal <- adrenal[,2]
  names(geneListadrenal) <- as.character(adrenal[,1])
  geneListadrenal <- sort(geneListadrenal, decreasing = TRUE)
  gene <- names(geneListadrenal)[abs(geneListadrenal) > 1]
  hsa05412  <- pathview(gene.data  = geneListadrenal, # hsa04151 PI3K
                       pathway.id = "hsa05412",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneListadrenal)), cpd=1),
                       out.suffix = "adrenal")
  
  
  
  
  