#' Create a heatmap for significant pathways
#'
#' Generates a heatmap based on the Over-Representation Analysis (ORA) results to visualize
#' significant pathways across communities. Each row in the heatmap represents a unique significant pathway,
#' while columns represent different communities. A value of 1 denotes significance, whereas 0 denotes non-significance.
#'
#' @param ora_results A list of ORA results, with each list element representing the results for a specific community.
#' @param significance_threshold A numeric threshold value for significance. Defaults to 0.05.
#' @importFrom pheatmap pheatmap
#' @return A heatmap visualizing the significant pathways across communities.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming `sample_ora_results` is a predefined list of ORA results.
#' ora_heatmap(ora_results = sample_ora_results, significance_threshold = 0.05)
#'}
ora_heatmap <- function(ora_results, significance_threshold = 0.05) {
  # 1. Extract significant pathways
  significant_gene_sets_list <- lapply(ora_results, function(result) {
    return(result$pathway[result$padj < significance_threshold])
  })

  # 2. Create binary matrix
  all_gene_sets <- unique(unlist(significant_gene_sets_list))
  binary_matrix <- matrix(0, nrow = length(all_gene_sets), ncol = length(ora_results),
                          dimnames = list(all_gene_sets, paste0("Community_", 1:length(ora_results))))

  for (i in 1:length(significant_gene_sets_list)) {
    binary_matrix[significant_gene_sets_list[[i]], i] <- 1
  }

  # 3. Plot the heatmap
  colors <- c("#377eb8", "#ff7f00")
  breaks <- c(0, 0.5, 1)

  pheatmap(binary_matrix,
           scale = "none",
           cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = T,
           color = colors,
           breaks = breaks,
           legend = T)
}
