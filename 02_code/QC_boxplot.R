
QC_boxplot <- function(data, data_group = NULL, value_colour = NULL, title) {

  library(ggplot2)
  library(tidyr)

  # Check if data_group is NULL or an empty data frame
  if (is.null(data_group) || (is.data.frame(data_group) && nrow(data_group) == 0)) {
    data_ggplot <- as.data.frame(t(data))
    data_ggplot$id <- rownames(data_ggplot)
    data_ggplot <- tidyr::gather(data_ggplot, key = "key", value = "value", -c("id"))

    ggplot(data_ggplot, aes(x = id, y = log2(value))) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "", title = title)
  } else {
    data_ggplot <- as.data.frame(t(data))
    data_ggplot$id <- rownames(data_ggplot)
    data_ggplot <- merge(data_group, data_ggplot, by = 'id')
    data_ggplot <- tidyr::gather(data_ggplot, key = "key", value = "value", -c("id", "group"))
    if (is.null(value_colour)) {
      stop("value_colour must be provided when data_group is not NULL")
    }

    if (length(value_colour) != length(unique(data_ggplot$group))) {
      warning("the length of value_colour is not equal to the length of data_group")
    }
    if (length(value_colour) < length(unique(data_ggplot$group))) {
      stop("the length of value_colour is less than the length of data_group")
    }

    data_ggplot <- data_ggplot[order(data_ggplot$group), ]
    data_ggplot$id <- factor(data_ggplot$id, levels = unique(data_ggplot$id))
    data_ggplot <- data_ggplot[order(data_ggplot$id), ]

    if (!is.factor(data_ggplot$group)) {
      data_ggplot$group <- factor(data_ggplot$group)
    }

    ggplot(data_ggplot, aes(x = id, y = log2(value), fill = group)) +
      geom_boxplot() +
      scale_fill_manual(values = value_colour) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "", title = title)
  }
}
