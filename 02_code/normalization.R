library(readr)
data_input <- read_delim("01_data/pg_120min.tsv", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
View(data_input)
data_input <- as.data.frame(data_input)
sample_id <- c("D1",
               "D2",
               "D3",
               "D4",
               "WT1",
               "WT2",
               "WT3",
               "WT4")
sample_group <- c(rep("D",4),
                  rep("WT",4))
value_colour <- c("WT" = "#4DBBD5FF",# control group
                  "D" = "#E64B35FF"# Experimental group
                  # "virginica" = "#00A087FF",# other group1
                  # "other group2" = "#3C5488FF"# other group2
                  )

data_group <- data.frame(id = sample_id,
                    group = sample_group)
rownames(data_input) <- data_input$Protein.Names
data_anno <- data_input[,!colnames(data_input) %in% sample_id]
data_matrix <- data_input[,colnames(data_input) %in% sample_id]
data_matrix[is.na(data_matrix)] <- 0
data_matrix <- data_matrix + 500
data_before <- data_matrix

# normalization -----------------------------------------------------------
data_after <- log2(data_matrix)
column_medians <- apply(data_after, 2, median, na.rm = TRUE)
target_median <- max(column_medians)  # 目标中位数
result_adjusted <- data_after

for (col in names(data_after)) {
  if (is.numeric(data_after[[col]])) {
    # 防止除以零的情况
    if (column_medians[col] != 0) {
      data_after[[col]] <- data_after[[col]] / column_medians[col] * target_median
    }
  }
}

column_medians_2 <- apply(data_after, 2, median, na.rm = TRUE)
data_after <- 2^data_after


# boxplot -----------------------------------------------------------------

QC_boxplot <- function(data, data_group, value_colour,title) {
  
  library(ggplot2)
  library(tidyr)
  
  data_ggplot <- as.data.frame(t(data))
  data_ggplot$id <- rownames(data_ggplot)
  data_ggplot <- merge(data_group,
                       data_ggplot,
                    by = 'id')
  data_ggplot <- tidyr::gather(data_ggplot,key = "key",
                               value = "value",
                               -c("id","group")
  )
  if (length(value_colour) != length(table(data_ggplot$group))){
    warning("the length of value_colour is not equal to the length of data_group")
  }
  if (length(value_colour) < length(table(data_ggplot$group))){
    stop("the length of value_colour is less than the length of data_group")
  }
  data_ggplot <- data_ggplot[order(data_ggplot$group),]
  data_ggplot$id <- factor(data_ggplot$id,levels = unique(data_ggplot$id))
  data_ggplot <- data_ggplot[order(data_ggplot$id),]
  if (!is.factor(data_ggplot$group)) {
    data_ggplot$group <- factor(data_ggplot$group)
  }
  ggplot(data_ggplot,aes(x = id,
                         y = log2(value + 500),
                         fill = group)
  ) +
    geom_boxplot() +
    scale_fill_manual(values = value_colour) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     colour = "black",
                                     size = 10),
          axis.text.y = element_text(hjust = 1,
                                     colour = "black",
                                     size = 10),
          plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "",
         title = title)

}
QC_boxplot(data_before,
           data_group = data_group,
           value_colour = value_colour,
           title = "rawdata")
QC_boxplot(data_after,
           data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")

# heatmap -----------------------------------------------------------------

QC_heatmap <- function(data,data_group,value_colour){
  library(tidyr)
  library(dplyr)
  data_pre <- as.data.frame(t(data))
  data_pre$id <- rownames(data_pre)
  data_pre <- merge(data_group,
                    data_pre,
                    by = 'id')
  rownames(data_pre) <- data_pre$id
   data_heatmap <- data_pre %>%
    subset(.,select = (-c(group,id))) %>%
    t() %>%
    as.data.frame() %>%
     {log2((. + 1))}
  library(pheatmap)
  # annotation_col requirements:
  # 1.the annotation_col must be a data frame
  # 2.the row names of annotation_col == the column names of data_heatmap
  # 3.the column names of annotation_col is the annotation legend name
  #
  # annotation_colors requirements:
  # 1.the annotation_colors must be a list.
  # 2.if the group is more than two, you can add by format.
  annotation_heatmap <- data_pre %>%
    select(group)
  
  annotation_colors <- list(group = value_colour)
  pheatmap(data_heatmap,
           show_rownames = F,
           annotation_col = annotation_heatmap,
           annotation_colors = annotation_colors,
           scale = "row")
}

QC_heatmap(data = data_before,
           data_group = data_group,
           value_colour = value_colour)

QC_heatmap(data = data_after,
           data_group = data_group,
           value_colour = value_colour)

# PCA ---------------------------------------------------------------------
QC_PCA <- function(data,data_group,value_colour){
  library(FactoMineR)
  library(factoextra)
  library(tidyr)
  library(dplyr)
  data_pre <- as.data.frame(t(data))
  data_pre$id <- rownames(data_pre)
  data_pre <- merge(data_group,
                    data_pre,
                    by = 'id')
  rownames(data_pre) <- data_pre$id
  group_list <- data_pre$group
  dat.pca <- data_pre %>%
    subset(.,select = -c(group,id)) %>%
    PCA(.,graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = c("text","point"), # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               palette = value_colour,
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
}

QC_PCA(data = data_before,
       data_group = data_group,
       value_colour = value_colour)

QC_PCA(data = data_after,
       data_group = data_group,
       value_colour = value_colour)
# DE ----------------------------------------------------------------------

library(limma)
group <- factor(data_group$group,levels = c("D","WT"))

limma_design<-model.matrix(~-1+group)

rownames(limma_design) <- data_group$id

limma_expr <- data

limma_expr <- log2(limma_expr)

contrast.matrix<-makeContrasts(contrasts = "groupD-groupWT",
                               levels = limma_design)  

fit <- lmFit(limma_expr,limma_design)   

fit1 <- contrasts.fit(fit, contrast.matrix)  

fit2 <- eBayes(fit1,0.01)    

tempOutput = topTable(fit2, 
                      coef = "groupD-groupWT", 
                      number = nrow(fit2),
                      lfc = log2(1),
                      adjust.method = "fdr")

dif <- tempOutput[tempOutput[,"adj.P.Val"] < 0.05,]
sd(tempOutput$logFC)

dif_logFC <- dif[abs(dif[,"logFC"]) > sd(tempOutput$logFC)*3,]

tempOutput$`Protein.Names` <- rownames(tempOutput) 
data_after$`Protein.Names` <- rownames(data_after)
result_merge <- merge(data_after,
                      tempOutput,
                      by = "Protein.Names")
result_merge <- merge(data_anno,
                      result_merge,
                      by = "Protein.Names")
tempOutput <- subset(tempOutput,select = -`Protein.Names`)
data_after <- subset(data_after,select = -`Protein.Names`)
write.csv(result_merge,file = "result.csv")
# volcano plot ------------------------------------------------------------

library(ggplot2)

library(ggrepel)

library(dplyr)

y <- rownames(dif_logFC)
anno_gene_symbol <- unlist(lapply(y,function(y) strsplit(as.character(y),"_")[[1]][1]))
Gene <- as.data.frame(anno_gene_symbol)
res_data <- tempOutput
data <- res_data[res_data[,"adj.P.Val"] <= 1,]

y <- rownames(data)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"_")[[1]][1]))
data$gene <- gene
# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
data$sig[data$`adj.P.Val` >= 0.05 | abs(data$logFC) < sd(tempOutput$logFC)*3] <- "Not"

data$sig[data$`adj.P.Val` < 0.05 & data$logFC >= sd(tempOutput$logFC)*3] <- "Up"

data$sig[data$`adj.P.Val` < 0.05 & data$logFC <= -sd(tempOutput$logFC)*3] <- "Down"
input <- data

volc <- ggplot(data = data, aes(x = logFC, 
                                y = -log10(`adj.P.Val`),
                                color = sig)) +
  theme_classic() +
  geom_point(alpha = 0.9) +
  theme_set(theme_bw()) + 
  theme(panel.grid = element_blank(),strip.text = element_blank()) +
  scale_color_manual(values = c("#4DBBD5","grey80","#E64B35")) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_vline(xintercept = c(-sd(tempOutput$logFC)*3,sd(tempOutput$logFC)*3),lty = 4,lwd = 0.6,alpha = 0.8) +
  labs(title = "D vs WT") +
  theme(plot.title = element_text(hjust = 0.5))


volc

volc + geom_text_repel(data = data[data$gene %in% Gene[,1],], aes(label=gene),size = 3,max.overlaps = 50)

# KEGG GO -----------------------------------------------------------------
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)

GeneSymbol <- subset(tempOutput,`adj.P.Val`< 0.05)
GeneSymbol <- subset(GeneSymbol,logFC > 0)

data <- data_anno[data_anno$Protein.Names %in% rownames(GeneSymbol) ,]
GeneSymbol <- data$Genes
gene.symbol.eg <- id2eg(ids = GeneSymbol, category = 'SYMBOL', org = 'Hs',na.rm = F) #如果是小鼠，就是Mm
gene.symbol.eg <- as.data.frame(gene.symbol.eg)

## GO分析 ----
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE) #如果是小鼠，就是Mm,out可选"ALL", "CC", "BP", "MF"
write.csv(kkd@result,file = "GO_up.csv",quote = F,row.names = F)
# 条形图
barplot(kkd,drop = T,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free',space = 'free') 
# 气泡图
dotplot(kkd,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free',space = 'free') 
#BiocManager::install('ggnewscale')
# 网络图
# cnetplot(kkd,categorySize = "star",foldChange = geneFC)

cnetplot(kkd,layout = 'star',foldChange = geneFC,colorEdge = T )


## KEGG分析 ----

kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID,organism = "hsa",keyType = "kegg",qvalueCutoff = 0.05)   #如果是小鼠，就是mouse
write.csv(kk@result,file = "KEGG_up.csv",quote = F,row.names = F)
kk@result$geneID <- gsub("/",",",kk@result$geneID)#更改基因间隔符号，否则没有logfc
# 条形图
barplot(kk,drop=T,showCategory = 10)
# 气泡图
dotplot(kk)
kk_plot <- as.data.frame(kk@result)
kk <- enrichKEGG(gene = mer$ENTREZID,organism = "hsa",keyType = "kegg",qvalueCutoff = 0.05) 
library(DOSE)
#如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
ego2<-setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# 网络图
cnetplot(ego2,layout = 'star',foldChange = geneFC,colorEdge = T)

