# DIA_report

此项目为DIANN下机数据处理

## 功能实现

此脚本通过读入DIA输出数据，完成intensity矫正，NA值填充，QC处理（boxplot，PCA plot和heatmap），差异分析和GO/KEGG通路富集分析。

提供示范示例**01_data/pg_120min.tsv**。

此项目QC功能实现来自于**Data-QC**项目。

运行脚本保存在`normalization.R`中，按步运行即可。

## function

### run_DE()

`run_DE()`函数使用如下：

``` r
# DE ----------------------------------------------------------------------
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
table(data_group$group)
# group 1为实验组
# group 2为对照组
group_1 <- "FL_9S"
group_2 <- "FL_non"
result_merge <- run_DE(data_fill = data_fill_normalization, # 填充后蛋白矩阵
                       data_group = data_group, # 分组信息，sampleid 为id列，group为group列
                       log2 = T, # 是否对数据进行log2运算
                       data_anno = data_anno, # 蛋白注释信息
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/")
```

`run_DE()`为**limma**运行流程。

运行后，该函数会自动在`dir`目录下生成`./group_vs_group2`目录，存放`result_merge.csv`结果文件。

同时该函数会自动将**data_fill**, **data_raw**, **data_anno**和**limma的差异结果**合并，方便后续检查差异结果。

### run_enrichment()

`run_enrichment()`函数使用如下：

``` r
# 根据上述limma DE结果，选择adj.P.Val或P.Value
GeneSymbol <- subset(result_merge,adj.P.Val< 0.05)
# GeneSymbol$Genes <- GeneSymbol$GN
y <- GeneSymbol$Genes
GeneSymbol$Genes <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))


result_enrich_pre <- run_enrichment_analysis(data = GeneSymbol, #column中Genes为gene symbol，logFC为差异倍数
                        OrgDb = "Hs", # human为Hs，mouse为Mm，Danio rerio为Dr
                        dir = paste0(DE_dir,
                                     group_1,"_vs_",group_2,"/"),
                        do_down_GO = F,   # 开启全部下调基因GO分析
                        do_down_KEGG = F, # 开启全部下调基因KEGG分析
                        do_up_GO = F,     # 开启全部下调基因GO分析
                        do_up_KEGG = F,   # 开启全部下调基因KEGG分析
                        do_all_GO = T,    # 开启全部基因GO分析
                        do_all_KEGG = T   # 开启全部基因KEGG分析
                                             )
```

`run_enrichment()`函数需存在Genes列为gene symbol，logFC为差异倍数。

`run_enrichment()`运行后会在`dir`目录下生成对应通路富集结果和log日志。如果存在某个通路富集无结果，请关闭对应通路富集分析开关后重新运行该函数。

### seqknn()

`seqknn()`函数使用如下：

``` R
# install.packages("multiUS")
library(multiUS)
# 输入原始蛋白表达矩阵
NA_ratio <-   colSums(is.na(data_input))/dim(data_input)[1]
NA_ratio <- as.data.frame(NA_ratio)
NA_ratio$Samples <- rownames(NA_ratio)

NA_ratio[NA_ratio$NA_ratio < 0.2,"group"] <- "Good"
NA_ratio[NA_ratio$NA_ratio < 0.5&NA_ratio$NA_ratio > 0.2,"group"] <- "OK"
NA_ratio[NA_ratio$NA_ratio < 0.8&NA_ratio$NA_ratio > 0.5,"group"] <- "Bad"
NA_ratio[NA_ratio$NA_ratio > 0.8,"group"] <- "Remove"
NA_ratio$group <- factor(NA_ratio$group,levels = c("Good","OK","Bad","Remove"))
library(ggplot2)
# 绘制每个样本的缺失值比例条形图
ggplot(NA_ratio,aes(x = NA_ratio, y = reorder(Samples,NA_ratio,decreasing = T), fill = group)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = sprintf("%.2f", NA_ratio)), hjust = -0.5) + 
  scale_fill_manual(values = c("Good" = "#62c882",
                               "OK" = "#b5e281",
                               "Bad" = "#febe80",
                               "Remove" = "#bb2022")) + 
  theme_bw() + 
  labs( x = "Number of missing rows",
        y = "Samples")
# 设定缺失值删除比例
cutoff_NA_ratio <- 0.5 #defult

NA_ratio_protein <- as.data.frame(rowSums(is.na(data_input))/dim(data_input)[2])
colnames(NA_ratio_protein) <- "NA_ratio_protein"
NA_ratio_protein$Protein.Group <- rownames(NA_ratio_protein)
NA_ratio_protein <- NA_ratio_protein[NA_ratio_protein$NA_ratio_protein < cutoff_NA_ratio,]
# NA_ratio_protein$protein_group <- rownames(NA_ratio_protein)
data <- data_input[rownames(data_input)%in%rownames(NA_ratio_protein),]
# 使用seqknn方法进行填充
data_fill <- multiUS::seqKNNimp(data = data,k = 10)
write.csv(data_fill,file = "01_rawdata/20241107_2ndDS/20241107_2ndDS/report.pg_matrix_fill_after.csv")
```
