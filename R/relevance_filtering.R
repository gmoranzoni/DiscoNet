#' Filter Nodes in an igraph Network Based on Relevance Score
#'
#' This function filters nodes in a given igraph network object based on the relevance score node attribute.
#' Nodes with relevance scores below the specified cutoff are removed from the network.
#'
#' @param network An igraph network object which must contain a vertex attribute named 'relevance_score'.
#' @param cutoff A numeric value specifying the threshold for the relevance score (default is 0).
#'        Nodes with a relevance score below this cutoff will be removed.
#' @return An igraph network object with nodes filtered based on the relevance score cutoff.
#'
#' @examples
#' \dontrun{
#' # Assuming 'net' is an igraph network object with a 'relevance_score' attribute for nodes:
#' filtered_net <- relevance_filtering(net, cutoff = 0.05)
#' }
#'
#' @export

relevance_filtering <- function(network, cutoff = 0) {
  # Ensure the input is an igraph object
  stopifnot(is.igraph(network))

  # Check if the 'relevance_score' attribute exists in the network nodes
  if (!"relevance_score" %in% vertex_attr_names(network)) {
    stop("Error: 'relevance_score' vertex attribute is missing.")
  }

  network_filtered <- delete_vertices(network, V(network)$relevance_score < cutoff)

  return(network_filtered)
}

