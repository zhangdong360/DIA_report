run_enrichment_analysis <- function(data, OrgDb, dir, 
                                   do_down_GO = TRUE, do_down_KEGG = TRUE,
                                   do_up_GO = TRUE, do_up_KEGG = TRUE,
                                   do_all_GO = FALSE, do_all_KEGG = FALSE) {
  library(clusterProfiler)
  library(pathview)
  library(ggplot2)
  
  # 加载对应的物种数据库
  if (OrgDb == "Hs") {
    library(org.Hs.eg.db)
  } else if (OrgDb == "Mm") {
    library(org.Mm.eg.db)
  } else if (OrgDb == "Dr") {
    library(org.Dr.eg.db)
  } else {
    stop("OrgDb参数必须是'Hs'(人类)、'Mm'(小鼠)或'Dr'(斑马鱼)")
  }
  
  # 初始进度提示
  message("\n[1/9] 开始富集分析流程...")
  message("=================================")
  message("分析设置:")
  message(paste(" - 物种:", switch(OrgDb, "Hs" = "人类", "Mm" = "小鼠", "Dr" = "斑马鱼")))
  message(paste(" - 下调基因GO分析:", ifelse(do_down_GO, "开启", "关闭")))
  message(paste(" - 下调基因KEGG分析:", ifelse(do_down_KEGG, "开启", "关闭")))
  message(paste(" - 上调基因GO分析:", ifelse(do_up_GO, "开启", "关闭")))
  message(paste(" - 上调基因KEGG分析:", ifelse(do_up_KEGG, "开启", "关闭")))
  message(paste(" - 全部基因GO分析:", ifelse(do_all_GO, "开启", "关闭")))
  message(paste(" - 全部基因KEGG分析:", ifelse(do_all_KEGG, "开启", "关闭")))
  message("=================================")
  
  # 检查输入数据是否有效
  message("[2/9] 检查输入数据有效性...")
  if (!all(c("Genes", "logFC") %in% colnames(data))) {
    stop("输入数据必须包含'Genes'和'logFC'列")
  }
  
  # 定义辅助函数
  function1 <- function(OrgDb, genesymbol, analysis_type, do_GO, do_KEGG) {
    if (length(genesymbol) == 0) {
      warning(paste(analysis_type, "基因列表为空"))
      return(list(enrichGO = NULL, enrichKEGG = NULL))
    }
    
    gene.symbol.eg <- NULL
    kkd <- NULL
    kk <- NULL
    
    tryCatch({
      # 基因ID转换
      message(paste0("  → 进行", analysis_type, "基因ID转换..."))
      if (OrgDb == "Hs") {
        gene.symbol.eg <- id2eg(ids = genesymbol, category = 'SYMBOL', org = 'Hs', na.rm = F)
      } else if (OrgDb == "Mm") {
        gene.symbol.eg <- id2eg(ids = genesymbol, category = 'SYMBOL', org = 'Mm', na.rm = F)
      } else if (OrgDb == "Dr") {
        gene.symbol.eg <- id2eg(ids = genesymbol, category = 'SYMBOL', org = 'dre', na.rm = F)
      }
      
      gene.symbol.eg <- as.data.frame(gene.symbol.eg)
      gene.symbol.eg <- na.omit(gene.symbol.eg)
      
      if (nrow(gene.symbol.eg) == 0) {
        warning(paste(analysis_type, "没有基因能成功转换为ENTREZID"))
        return(list(enrichGO = NULL, enrichKEGG = NULL))
      }
      
      message(paste0("  √ ", analysis_type, "基因ID转换完成 (", nrow(gene.symbol.eg), "个基因转换成功)"))
      
      # GO分析
      if (do_GO) {
        message(paste0("  → 进行", analysis_type, "GO富集分析..."))
        if (OrgDb == "Hs") {
          kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                         OrgDb = "org.Hs.eg.db",
                         ont = "ALL",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
        } else if (OrgDb == "Mm") {
          kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                         OrgDb = "org.Mm.eg.db",
                         ont = "ALL",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
        } else if (OrgDb == "Dr") {
          kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                         OrgDb = "org.Dr.eg.db",
                         ont = "ALL",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
        }
        
        if (!is.null(kkd) && nrow(kkd@result) == 0) {
          warning(paste(analysis_type, "GO富集分析没有显著结果"))
          kkd <- NULL
        } else if (!is.null(kkd)) {
          message(paste0("  √ ", analysis_type, "GO富集分析完成 (", nrow(kkd@result), "条显著通路)"))
        }
      } else {
        message(paste0("  × 跳过", analysis_type, "GO富集分析(用户设置)"))
      }
      
      # KEGG分析
      if (do_KEGG) {
        message(paste0("  → 进行", analysis_type, "KEGG富集分析..."))
        if (OrgDb == "Hs") {
          kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", keyType = "kegg", qvalueCutoff = 0.05)
        } else if (OrgDb == "Mm") {
          kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "mmu", keyType = "kegg", qvalueCutoff = 0.05)
        } else if (OrgDb == "Dr") {
          kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "dre", keyType = "kegg", qvalueCutoff = 0.05)
        }
        
        if (!is.null(kk) && nrow(kk@result) == 0) {
          warning(paste(analysis_type, "KEGG富集分析没有显著结果"))
          kk <- NULL
        } else if (!is.null(kk)) {
          message(paste0("  √ ", analysis_type, "KEGG富集分析完成 (", nrow(kk@result), "条显著通路)"))
        }
      } else {
        message(paste0("  × 跳过", analysis_type, "KEGG富集分析(用户设置)"))
      }
      
      result <- list(enrichGO = kkd, enrichKEGG = kk)
      return(result)
    }, error = function(e) {
      warning(paste(analysis_type, "富集分析过程中出错:", e$message))
      return(list(enrichGO = NULL, enrichKEGG = NULL))
    })
  }
  
  # 创建目录
  message("[3/9] 创建输出目录...")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(paste("√ 目录创建成功:", dir))
  } else {
    message(paste("√ 目录已存在:", dir))
  }
  
  # down-regulated genes
  message("\n[4/9] 分析下调基因...")
  message("-----------------------------")
  down_genes <- subset(data, logFC < 0)
  if (nrow(down_genes) == 0) {
    warning("没有下调基因，跳过下调基因分析")
  } else {
    message(paste("↓ 发现", nrow(down_genes), "个下调基因"))
    
    # 只有当至少一个下调分析开启时才执行
    if (do_down_GO || do_down_KEGG) {
      result_down <- function1(OrgDb, down_genes$Genes, "下调基因", do_down_GO, do_down_KEGG)
      kkd_down <- result_down$enrichGO
      kk_down <- result_down$enrichKEGG
      
      # 处理下调基因GO结果
      if (do_down_GO && !is.null(kkd_down)) {
        tryCatch({
          message("  → 输出下调基因GO分析结果...")
          write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
          p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
          print(p1)
          dev.off()
          message("  √ 下调基因GO结果输出完成")
        }, error = function(e) {
          warning(paste("下调基因GO分析结果输出失败:", e$message))
        })
      }
      
      # 处理下调基因KEGG结果
      if (do_down_KEGG && !is.null(kk_down)) {
        tryCatch({
          message("  → 输出下调基因KEGG分析结果...")
          write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
          p2 <- dotplot(kk_down)
          print(p2)
          dev.off()
          message("  √ 下调基因KEGG结果输出完成")
        }, error = function(e) {
          warning(paste("下调基因KEGG分析结果输出失败:", e$message))
        })
      }
    } else {
      message("  × 下调基因分析已全部关闭(用户设置)")
    }
  }
  
  # up-regulated genes
  message("\n[5/9] 分析上调基因...")
  message("-----------------------------")
  up_genes <- subset(data, logFC > 0)
  if (nrow(up_genes) == 0) {
    warning("没有上调基因，跳过上调基因分析")
  } else {
    message(paste("↑ 发现", nrow(up_genes), "个上调基因"))
    
    # 只有当至少一个上调分析开启时才执行
    if (do_up_GO || do_up_KEGG) {
      result_up <- function1(OrgDb, up_genes$Genes, "上调基因", do_up_GO, do_up_KEGG)
      kkd_up <- result_up$enrichGO
      kk_up <- result_up$enrichKEGG
      
      # 处理上调基因GO结果
      if (do_up_GO && !is.null(kkd_up)) {
        tryCatch({
          message("  → 输出上调基因GO分析结果...")
          write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
          p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
          print(p3)
          dev.off()
          message("  √ 上调基因GO结果输出完成")
        }, error = function(e) {
          warning(paste("上调基因GO分析结果输出失败:", e$message))
        })
      }
      
      # 处理上调基因KEGG结果
      if (do_up_KEGG && !is.null(kk_up)) {
        tryCatch({
          message("  → 输出上调基因KEGG分析结果...")
          write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
          p4 <- dotplot(kk_up)
          print(p4)
          dev.off()
          message("  √ 上调基因KEGG结果输出完成")
        }, error = function(e) {
          warning(paste("上调基因KEGG分析结果输出失败:", e$message))
        })
      }
    } else {
      message("  × 上调基因分析已全部关闭(用户设置)")
    }
  }
  
  # all genes analysis
  message("\n[6/9] 分析全部基因...")
  message("-----------------------------")
  all_genes <- data$Genes
  if (length(all_genes) == 0) {
    warning("没有基因可用于分析")
  } else {
    message(paste("发现", length(all_genes), "个基因"))
    
    # 只有当至少一个全部基因分析开启时才执行
    if (do_all_GO || do_all_KEGG) {
      result_all <- function1(OrgDb, all_genes, "全部基因", do_all_GO, do_all_KEGG)
      kkd_all <- result_all$enrichGO
      kk_all <- result_all$enrichKEGG
      
      # 处理全部基因GO结果
      if (do_all_GO && !is.null(kkd_all)) {
        tryCatch({
          message("  → 输出全部基因GO分析结果...")
          write.csv(kkd_all@result, file = paste0(dir, "/GO_all.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/GO_all.pdf"), width = 6, height = 7)
          p5 <- dotplot(kkd_all, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
          print(p5)
          dev.off()
          message("  √ 全部基因GO结果输出完成")
        }, error = function(e) {
          warning(paste("全部基因GO分析结果输出失败:", e$message))
        })
      }
      
      # 处理全部基因KEGG结果
      if (do_all_KEGG && !is.null(kk_all)) {
        tryCatch({
          message("  → 输出全部基因KEGG分析结果...")
          write.csv(kk_all@result, file = paste0(dir, "/KEGG_all.csv"), quote = F, row.names = F)
          pdf(file = paste0(dir, "/KEGG_all.pdf"), width = 6, height = 5)
          p6 <- dotplot(kk_all)
          print(p6)
          dev.off()
          message("  √ 全部基因KEGG结果输出完成")
        }, error = function(e) {
          warning(paste("全部基因KEGG分析结果输出失败:", e$message))
        })
      }
    } else {
      message("  × 全部基因分析已全部关闭(用户设置)")
    }
  }
  
  # 结果汇总
  message("\n[7/9] 分析结果汇总:")
  message("-----------------------------")
  if (exists("kkd_down") && !is.null(kkd_down)) message(paste("下调基因GO通路:", nrow(kkd_down@result), "条"))
  if (exists("kk_down") && !is.null(kk_down)) message(paste("下调基因KEGG通路:", nrow(kk_down@result), "条"))
  if (exists("kkd_up") && !is.null(kkd_up)) message(paste("上调基因GO通路:", nrow(kkd_up@result), "条"))
  if (exists("kk_up") && !is.null(kk_up)) message(paste("上调基因KEGG通路:", nrow(kk_up@result), "条"))
  if (exists("kkd_all") && !is.null(kkd_all)) message(paste("全部基因GO通路:", nrow(kkd_all@result), "条"))
  if (exists("kk_all") && !is.null(kk_all)) message(paste("全部基因KEGG通路:", nrow(kk_all@result), "条"))
  
  # 保存分析日志
  message("[8/9] 保存分析日志...")
  sink(file = paste0(dir, "/analysis_log.txt"))
  cat("富集分析日志\n")
  cat("=============\n")
  cat("分析时间:", date(), "\n\n")
  cat("分析设置:\n")
  cat(paste(" - 物种:", switch(OrgDb, "Hs" = "人类", "Mm" = "小鼠", "Dr" = "斑马鱼"), "\n"))
  cat(paste(" - 下调基因GO分析:", ifelse(do_down_GO, "开启", "关闭"), "\n"))
  cat(paste(" - 下调基因KEGG分析:", ifelse(do_down_KEGG, "开启", "关闭"), "\n"))
  cat(paste(" - 上调基因GO分析:", ifelse(do_up_GO, "开启", "关闭"), "\n"))
  cat(paste(" - 上调基因KEGG分析:", ifelse(do_up_KEGG, "开启", "关闭"), "\n"))
  cat(paste(" - 全部基因GO分析:", ifelse(do_all_GO, "开启", "关闭"), "\n"))
  cat(paste(" - 全部基因KEGG分析:", ifelse(do_all_KEGG, "开启", "关闭"), "\n\n"))
  
  cat("基因统计:\n")
  cat(paste(" - 下调基因数量:", ifelse(exists("down_genes"), nrow(down_genes), 0), "\n"))
  cat(paste(" - 上调基因数量:", ifelse(exists("up_genes"), nrow(up_genes), 0), "\n"))
  cat(paste(" - 全部基因数量:", length(all_genes), "\n\n"))
  
  cat("富集结果汇总:\n")
  if (exists("kkd_down") && !is.null(kkd_down)) cat(paste(" - 下调基因GO通路:", nrow(kkd_down@result), "\n"))
  if (exists("kk_down") && !is.null(kk_down)) cat(paste(" - 下调基因KEGG通路:", nrow(kk_down@result), "\n"))
  if (exists("kkd_up") && !is.null(kkd_up)) cat(paste(" - 上调基因GO通路:", nrow(kkd_up@result), "\n"))
  if (exists("kk_up") && !is.null(kk_up)) cat(paste(" - 上调基因KEGG通路:", nrow(kk_up@result), "\n"))
  if (exists("kkd_all") && !is.null(kkd_all)) cat(paste(" - 全部基因GO通路:", nrow(kkd_all@result), "\n"))
  if (exists("kk_all") && !is.null(kk_all)) cat(paste(" - 全部基因KEGG通路:", nrow(kk_all@result), "\n"))
  sink()
  message("√ 分析日志保存完成")
  
  # 完成提示
  message("\n[9/9] 富集分析完成!")
  message("=================================")
  message(paste("结果已保存到:", normalizePath(dir)))
  message("请检查警告信息(如果有)\n")
  
  # 返回结果(可选)
  result_list <- list()
  if (exists("kkd_down")) result_list$down_GO <- kkd_down
  if (exists("kk_down")) result_list$down_KEGG <- kk_down
  if (exists("kkd_up")) result_list$up_GO <- kkd_up
  if (exists("kk_up")) result_list$up_KEGG <- kk_up
  if (exists("kkd_all")) result_list$all_GO <- kkd_all
  if (exists("kk_all")) result_list$all_KEGG <- kk_all
  
  return(result_list)
}