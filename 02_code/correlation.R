# 读取数据
library(readr)
all_NA.pg <- read_csv("01_data/2501all/all_NA.pg.csv")
all_NA.pg <- as.data.frame(all_NA.pg)
rownames(all_NA.pg) <- all_NA.pg$...1
all_NA.pg <- subset(all_NA.pg, select = -c(`...1`))

# 计算相关性矩阵
result_cor <- cor(all_NA.pg, use = "complete.obs", method = "pearson")

# 检查相关性矩阵是否计算成功
print(result_cor)

# 绘制相关性图
library(corrplot)
corrplot::corrplot.mixed(result_cor,tl.col = "black",is.corr = F,tl.pos = "lt",
                         # order = "hclust",hclust.method = "complete",
                         upper = "number",lower = "circle")
