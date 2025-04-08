# data input ----
## function ----
source("02_code/QC_PCA.R")
source("02_code/QC_boxplot.R")
source("02_code/QC_heatmap.R")
source("02_code/run_DE.R")
source("02_code/run_enrichment_analysis.R")
## group input ----
library(readxl)
library(ggplot2)
# 导入分组信息
data_group <- read_excel("./01_data/group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
# 配色设置
value_colour <- c("D" = "#E64B35FF",# Experimental group
                  "WT" = "#4DBBD5FF"# other group1
)
rownames(data_group) <- data_group$id
## DIA matrix input ----
library(readr)
library(dplyr)
data_input <- read_delim("./01_data/pg_120min.tsv",
                         delim = "\t", escape_double = FALSE,
                         trim_ws = TRUE)
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Protein.Group
# 保留前五列注释
data_anno <- data_input[,1:5]
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
data_input <- data_input[,-1:-5]

## NAguideR ----
# 是否进行log2运算？
data_input[,-1:-5] <- log2(data_input[,-1:-5])
# 是否进行median normalization运算？
median_values <- apply(data_input[,-1:-5], 2, median, na.rm = TRUE)
data_input[,-1:-5] <- data_input[,-1:-5] %>%#除以每列中位数
  mutate(across(everything() , ~ ./median(., na.rm = TRUE)))
# 使用此输出进行NA填充，填充网站 https://www.omicsolution.org/wukong/NAguideR/#
write.csv(data_input,file = "01_rawdata/report.pg_matrix_fill_before.csv")

# 将填充后的数据导入
data_fill <- read_csv("01_rawdata/report.pg_matrix_fill_after.csv")
data_fill <- as.data.frame(data_fill)
rownames(data_fill) <- data_fill$Protein.Group
data_fill <- subset(data_fill,select = -c(`Protein.Group`))

# 如果进行了median normalization运算，运行下面的函数
data_fill <- sweep(data_fill, 2, median_values, `*`)
# 如果进行log2处理，运行下面的函数
data_fill <- 2 ^ data_fill
## 10minimum fill ----
intensity_all_protein <- data.frame(Protein.Group = data_input$Protein.Group,
                        intensity = rowMeans(data_input[,-1:-5],na.rm = T))
intensity = rowMeans(data_input[,-1:-5],na.rm = T)
intensity <- na.omit(intensity)
intensity <- sort(intensity)
intensity_10 <- quantile(intensity, 0.1)
intensity_10
log2(intensity_10)
ggplot(intensity, aes(x = log2(intensity))) + 
  geom_density(color = "black", fill = "gray")

# 将第6列到最后一列的所有NA值替换为0
data_input[,-1:-5][is.na(data_input[,-1:-5])] <- 0

# 将第6列到最后一列的所有值加100000
data_input[,-1:-5] <- data_input[,-1:-5] + 100000
data_fill <- data_input[,-1:-5]
# normalization -----------------------------------------------------------
# 设置输出目录
dir <- "03_result/"
## intensity normalization ----
data_before <- log2(data_fill)
# 计算校正前各样本的intensity median
column_medians <- apply(data_before, 2, median, na.rm = TRUE)
column_medians
# 目标中位数
target_median <- max(column_medians)
data_after <- data_before
for (col in names(data_after)) {
  if (is.numeric(data_after[[col]])) {
    # 防止除以零的情况
    if (column_medians[col] != 0) {
      data_after[[col]] <- data_after[[col]] / column_medians[col] * target_median
    }
  }
}
# 验证校正后各样本intensity median是否一致
column_medians_2 <- apply(data_after, 2, median, na.rm = TRUE)
column_medians_2
# 返回log2之前的数据
data_fill_normalization <- 2 ^ data_after
rm(data_before)

write.csv(data_fill_normalization,file = "01_rawdata/report.pg_matrix_fill_normalization.csv")


# boxplot -----------------------------------------------------------------
pdf(file = paste0(dir,"QC_boxplot_before.pdf"),
    width = 6,
    height = 4)

QC_boxplot(data_before,data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_fill_normalization,data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
# heatmap -----------------------------------------------------------------

pdf(file = paste0(dir,"QC_heatmap_before.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_before,data_group = data_group,
           value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_fill_normalization,data_group = data_group,
           value_colour = value_colour)
dev.off()

# PCA ---------------------------------------------------------------------
pdf(file = paste0(dir,"QC_pca_before.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = data_before,
       data_group = data_group,
       value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = data_fill_normalization,
       data_group = data_group,
       value_colour = value_colour)
dev.off()
# DE ----------------------------------------------------------------------
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
table(data_group$group)
# group 1为实验组
# group 2为对照组
group_1 <- "FL_9S"
group_2 <- "FL_non"
result_merge <- run_DE(data_fill = data_fill_normalization,
                       data_group = data_group,
                       log2 = T,
                       data_anno = data_anno,
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/")

# volcano plot ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

dif_logFC_up <- subset(result_merge,result_merge$logFC > 0)
dif_logFC_up <- dif_logFC_up[order(dif_logFC_up$logFC,decreasing = T),]
dif_logFC_down <- subset(result_merge,result_merge$logFC < 0)
dif_logFC_down <- dif_logFC_down[order(dif_logFC_down$logFC),]

# y <- c(rownames(dif_logFC_up)[1:10],rownames(dif_logFC_down)[1:10])
# y <- na.omit(y)

res_data <- result_merge
data <- res_data[res_data[,"P.Value"] <= 1,]

y <- result_merge$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
data$gene <- gene
# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
data$sig[data$P.Value >= 0.05 | abs(data$logFC) < 0] <- "Not"

data$sig[data$P.Value < 0.05 & data$logFC >= 0] <- "Up"

data$sig[data$P.Value < 0.05 & data$logFC <= -0] <- "Down"
input <- data
library(ggrepel)
library(ggplot2)
volc <- ggplot(data = data, aes(x = logFC,
                                y = -log10(P.Value),
                                color = sig)) +
  geom_point(alpha = 0.9) +  theme_classic() +
  theme(panel.grid = element_blank(),strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#4DBBD5","grey80","#E64B35")) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +

  labs(title = paste0(group_1,"-",group_2))
volc



## 标记火山图 ----
library(readxl)
list <- read_xlsx("01_rawdata/A list of positive control proteins.xlsx")
y <- list$`UniProt accession`
anno_gene_symbol <- unlist(lapply(y,function(y) strsplit(as.character(y),"_")[[1]][1]))
Gene <- as.data.frame(anno_gene_symbol)
library(readr)
rownames(result_merge) <- result_merge$Protein.Names

volc +
  geom_text_repel(data = data[data$Protein.Group %in% Gene[,1],],
                       aes(label=gene),
                  size = 3,
                  color = "black",
                  max.overlaps = 50) +
  geom_point(data = data[data$Protein.Group %in% Gene[,1],],
             alpha = 0.9,color = "black")
ggsave(plot = last_plot(),filename = paste0(dir,"volc.pdf"),
       height = 5,
       width = 6)

# KEGG GO -----------------------------------------------------------------
GeneSymbol <- subset(result_merge,adj.P.Val< 0.05)
# GeneSymbol$Genes <- GeneSymbol$GN
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))



# 样本为人，OrgDb为Hs
# 样本为小鼠，OrgDb为Mm
run_enrichment_analysis(data = GeneSymbol,
                        OrgDb = "Hs",
                        dir = paste0("03_result/DE/",
                                     group_1,"_vs_",group_2,"/"))

