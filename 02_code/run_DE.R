run_DE <- function(data, data_group, data_anno=NULL,log2=TRUE,group_1,group_2,dir = getwd()) {
  library(limma)
  # check data
  # 设置输出目录，并生成目录
  # 检查行名是否一致
  if (!is.null(data_anno)){
    if (!all(rownames(data) == rownames(data_anno))) {
      stop("Error: The row names of 'data' and 'data_anno' do not match.")
    }
  }
  output_dir <- paste0(dir, "/", group_1, "_vs_", group_2)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # 筛选对应分组
  limma_group <- subset(data_group,data_group$group %in% c(group_1,group_2))
  group <- factor(limma_group$group,levels = c(group_1,group_2))# "D","WT"

  limma_design <- model.matrix(~-1+group)

  rownames(limma_design) <- limma_group$id

  limma_expr <- data[,colnames(data) %in% limma_group$id]
  if(log2 == TRUE){
    limma_expr <- log2(limma_expr)
  }

  contrast.matrix <- makeContrasts(contrasts = paste0("group",levels(group)[1],"-","group",levels(group)[2]),
                                   levels = limma_design)

  fit <- lmFit(limma_expr,limma_design)

  fit1 <- contrasts.fit(fit, contrast.matrix)

  fit2 <- eBayes(fit1,0.01)

  tempOutput = topTable(fit2,
                        coef = paste0("group",levels(group)[1],"-","group",levels(group)[2]),
                        number = nrow(fit2),
                        lfc = log2(1),
                        adjust.method = "fdr")
  result_merge <- merge(data,
                        tempOutput,
                        by = 0)
  rownames(result_merge) <- result_merge$Row.names
  if (!is.null(data_anno)){
    result_merge <- subset(result_merge,select = -c(`Row.names`))
    result_merge <- merge(data_anno,
                          result_merge,
                          by = 0)
  }
  result_merge <- subset(result_merge,select = -c(`Row.names`))
  write.csv(result_merge,file = paste0(output_dir,"/result_DE.csv"))
  return(result_merge)
}
