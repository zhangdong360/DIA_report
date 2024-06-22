QC_PCA <- function(data, data_group = NULL, value_colour = NULL) {
  library(FactoMineR)
  library(factoextra)
  library(tidyr)
  library(dplyr)

  # Check if data_group is NULL or an empty data frame
  if (is.null(data_group) || (is.data.frame(data_group) && nrow(data_group) == 0)) {
    data_pre <- as.data.frame(scale(t(data)))
    dat.pca <- PCA(data_pre, graph = FALSE)
    fviz_pca_ind(dat.pca,
                 geom.ind = c("text", "point"), # show points only (but not "text")
                 addEllipses = FALSE # No concentration ellipses
    )
  } else {
    data_pre <- as.data.frame(t(data))
    data_pre$id <- rownames(data_pre)
    data_pre <- merge(data_group, data_pre, by = 'id')
    rownames(data_pre) <- data_pre$id
    group_list <- data_pre$group

    # Check if value_colour is provided
    if (is.null(value_colour)) {
      stop("value_colour must be provided when data_group is not NULL")
    }

    dat.pca <- data_pre %>%
      dplyr::select(-c(group, id)) %>%
      PCA(graph = FALSE)
    fviz_pca_ind(dat.pca,
                 geom.ind = c("text", "point"), # show points only (but not "text")
                 col.ind = group_list, # color by groups
                 palette = value_colour,
                 addEllipses = TRUE, # Concentration ellipses
                 legend.title = "Groups"
    ) +
      theme_classic()
  }
}
