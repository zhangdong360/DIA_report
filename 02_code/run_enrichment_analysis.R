run_enrichment_analysis <- function(data, OrgDb, dir) {
  library(clusterProfiler)
  library(pathview)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)

  # 定义辅助函数
  function1 <- function(OrgDb, genesymbol) {
    if (OrgDb == "Hs") {
      print("OrgDb: Human")
      gene.symbol.eg <- id2eg(ids = genesymbol, category = 'SYMBOL', org = 'Hs', na.rm = F)
    } else if (OrgDb == "Mm") {
      print("OrgDb: Mouse")
      gene.symbol.eg <- id2eg(ids = genesymbol, category = 'SYMBOL', org = 'Mm', na.rm = F)
    }
    gene.symbol.eg <- as.data.frame(gene.symbol.eg)
    gene.symbol.eg <- na.omit(gene.symbol.eg)

    # GO分析
    if (OrgDb == "Hs") {
      kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                      OrgDb = "org.Hs.eg.db",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
      kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", keyType = "kegg", qvalueCutoff = 0.05)
    } else if (OrgDb == "Mm") {
      kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                      OrgDb = "org.Mm.eg.db",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
      kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "mmu", keyType = "kegg", qvalueCutoff = 0.05)
    }

    result <- list(enrichGO = kkd, enrichKEGG = kk)
    return(result)
  }

  # 创建目录
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # down-regulated genes
  down_genes <- subset(data, logFC < 0)
  result_down <- function1(OrgDb, down_genes$Genes)
  kkd_down <- result_down$enrichGO
  kk_down <- result_down$enrichKEGG

  write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
  p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
  print(p1)
  dev.off()

  write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
  p2 <- dotplot(kk_down)
  print(p2)
  dev.off()

  # up-regulated genes
  up_genes <- subset(data, logFC > 0)
  result_up <- function1(OrgDb, up_genes$Genes)
  kkd_up <- result_up$enrichGO
  kk_up <- result_up$enrichKEGG

  write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
  p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
  print(p3)
  dev.off()

  write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
  p4 <- dotplot(kk_up)
  print(p4)
  dev.off()
}
