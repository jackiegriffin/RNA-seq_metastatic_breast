# libraries ----
  library(clusterProfiler)
  library(magrittr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)

# bone data ----
  d_bone <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/bone.2.11.2021.csv")

  # Format data ----
  geneListbone <- d_bone[,2]
  names(geneListbone) <- as.character(d_bone[,1])## feature 2: named vector
  geneListbone <- sort(geneListbone, decreasing = TRUE)## feature 3: decreasing order
  genebone <- names(geneListbone)[abs(geneListbone) > 2]
  gene.df.bone <- bitr(genebone, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                   OrgDb = org.Hs.eg.db)
  head(gene.df.bone)

# KEGG ----
  
  
  enrich_KEGG_bone <- enrichKEGG(gene = genebone, organism='hsa',pvalueCutoff = 0.05)
  KEGG_bone <- dotplot(enrich_KEGG_bone, showCategory=20) + ggtitle("KEGG ORA (Bone)")
  plot(KEGG_bone)
  # plot_grid(KEGG_bone, KEGG_brain, cols = 1)
  KEGG_CNET_BONE <- cnetplot(enrich_KEGG_bone, foldChange=geneListbone)
  plot(KEGG_CNET_BONE)
# GO ORA ---
ego_bone <- enrichGO(gene          = genebone,
                universe      = names(geneListbone),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                # qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_bone)

GO_ORA_bar_bone <- barplot(ego_bone, showCategory=40)
plot(GO_ORA_bar_bone)
GO_ORA_cnet_bone <- cnetplot(ego_bone, foldChange=geneListbone)
plot(GO_ORA_cnet_bone)



# GO gene set enrighment ----
ggo_bone <- groupGO(gene     = genebone,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo_bone)

go_group_bone <-barplot(ggo_bone, showCategory=40)
plot(go_group_bone)




# ------------------
ego3_bone <- gseGO(geneList     = geneListbone,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              # pvalueCutoff = 0.05,
              verbose      = FALSE)
head(ego3_bone)


ego3_bone_bone <- dotplot(ego3_bone, showCategory=20) + ggtitle("gseGO (Bone)")
plot(ego3_bone_bone)
