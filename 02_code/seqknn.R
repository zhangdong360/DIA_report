# install.packages("multiUS")
library(multiUS)
NA_ratio <-   colSums(is.na(data_input))/dim(data_input)[1]
NA_ratio <- as.data.frame(NA_ratio)
NA_ratio$Samples <- rownames(NA_ratio)

NA_ratio[NA_ratio$NA_ratio < 0.2,"group"] <- "Good"
NA_ratio[NA_ratio$NA_ratio < 0.5&NA_ratio$NA_ratio > 0.2,"group"] <- "OK"
NA_ratio[NA_ratio$NA_ratio < 0.8&NA_ratio$NA_ratio > 0.5,"group"] <- "Bad"
NA_ratio[NA_ratio$NA_ratio > 0.8,"group"] <- "Remove"
NA_ratio$group <- factor(NA_ratio$group,levels = c("Good","OK","Bad","Remove"))
library(ggplot2)
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

cutoff_NA_ratio <- 0.5 #defult

NA_ratio_protein <- as.data.frame(rowSums(is.na(data_input))/dim(data_input)[2])
colnames(NA_ratio_protein) <- "NA_ratio_protein"
NA_ratio_protein$Protein.Group <- rownames(NA_ratio_protein)
NA_ratio_protein <- NA_ratio_protein[NA_ratio_protein$NA_ratio_protein < cutoff_NA_ratio,]
# NA_ratio_protein$protein_group <- rownames(NA_ratio_protein)
data <- data_input[rownames(data_input)%in%rownames(NA_ratio_protein),]
data_fill <- multiUS::seqKNNimp(data = data,k = 10)
write.csv(data_fill,file = "01_rawdata/20241107_2ndDS/20241107_2ndDS/report.pg_matrix_fill_after.csv")
