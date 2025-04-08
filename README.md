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

```R
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

```R
result_enrich_pre <- run_enrichment_analysis(data = GeneSymbol, #Genes column为gene symbol，logFC为差异倍数
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
