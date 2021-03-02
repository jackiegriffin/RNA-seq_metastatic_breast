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

# Format data ----
geneListadrenal <- d_adrenal[,2]
names(geneListadrenal) <- as.character(d_adrenal[,1])## feature 2: named vector
geneListadrenal <- sort(geneListadrenal, decreasing = TRUE)## feature 3: decreasing order
geneadrenal <- names(geneListadrenal)[abs(geneListadrenal) > 1]
gene.df.adrenal <- bitr(geneadrenal, fromType = "ENTREZID",
                     toType = c("ENSEMBL", "SYMBOL"),
                     OrgDb = org.Hs.eg.db)
head(gene.df.adrenal)

# KEGG ----


enrich_KEGG_adrenal <- enrichKEGG(gene = geneadrenal, organism='hsa',pvalueCutoff = 0.05)
KEGG_adrenal <- dotplot(enrich_KEGG_adrenal, showCategory=20) + ggtitle("KEGG ORA for adrenal metastases")
plot(KEGG_adrenal)
plot_grid(KEGG_adrenal, KEGG_brain, KEGG_bone, ncol = 1)



ego3_adrenal <- gseGO(geneList     = geneListadrenal,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "CC",
                    nPerm        = 1000,
                    minGSSize    = 100,
                    maxGSSize    = 500,
                    # pvalueCutoff = 0.05,
                    verbose      = FALSE)
head(ego3_adrenal)


ego3_adrenal_dot <- dotplot(ego3_adrenal, showCategory=20) + ggtitle("gseGO (Adrenal)")
plot(ego3_adrenal_dot)

plot_grid(ego3_adrenal_dot, ego3_brain_dot, ego3_bone_bone, ncol = 2)













# GO ORA ---
ego_adrenal <- enrichGO(gene          = gene,
                     universe      = names(geneListadrenal),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(ego_adrenal)

GO_ORA_bar_adrenal <- barplot(ego_adrenal, showCategory=40)
plot(GO_ORA_bar_adrenal)
GO_ORA_cnet_adrenal <- cnetplot(ego_adrenal, foldChange=geneListadrenal)
plot(GO_ORA_cnet_adrenal)



# GO gene set enrighment ----
ggo_adrenal <- groupGO(gene     = gene,
                    OrgDb    = org.Hs.eg.db,
                    ont      = "CC",
                    level    = 3,
                    readable = TRUE)

head(ggo_adrenal)

go_group_adrenal <-barplot(ggo_adrenal, showCategory=40)
plot(go_group_adrenal)




