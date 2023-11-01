#' Convert Network Data to Igraph Object with Attributes
#'
#' This function creates an igraph object from the provided network data. It assigns
#' the attributes stored in `node_table` to the vertices (vertex attributes) of the igraph object.
#' The function also offers the possibility to set a relevance score cutoff to retrieve a network
#' containing only genes/proteins with a relevance score higher than the specified cutoff.
#' In the absence of a `node_table` or `edge_table`, the input network is converted into an igraph
#' object with no attributes.
#'
#' @param network A dataframe representing the network, expected to have columns that can be used as edges.
#' @param node_table A dataframe or NULL, representing node attributes. Defaults to NULL.
#' @param edge_table A dataframe or NULL, representing edge attributes. Defaults to NULL.
#' @param relevance_score_cutoff A numeric value representing the minimum relevance score for nodes to be included in the output graph. Defaults to 0.
#'
#' @return An igraph object representing the network, with specified attributes.
#' @import igraph
#' @export
#'
#' @examples
#' \dontrun{
#' # Using example data frames for network, node_table, and edge_table (to be provided for practical examples).
#'
#' network_df <- data.frame(node1 = c('A', 'B', 'C'), node2 = c('B', 'C', 'D'))
#' node_attributes <- data.frame(nodes = c('A', 'B', 'C', 'D'), attribute1 = c(1, 2, 3, 4), attribute2 = c(5, 6, 7, 8))
#' edge_attributes <- data.frame(edge_attr = c('X', 'Y', 'Z')) # Dummy example, adapt accordingly
#'
#' graph <- network_igraph(network = network_df, node_table = node_attributes, edge_table = edge_attributes, relevance_score_cutoff = 2)
#'}
network_igraph <- function(network, node_table = NULL, edge_table = NULL, relevance_score_cutoff = 0) {

  network_graph <- graph_from_data_frame(d = network, directed = F)

  ## Only if the network is of order 1
  # Assign the vertex attributes
  if("ecc_rel_score" %in% colnames(network) & "seed" %in% colnames(network)) {
  V(network_graph)$ecc_rel_score <- network$ecc_rel_score[match(V(network_graph)$name, network$nodes)]
  V(network_graph)$seed <- network$seed[match(V(network_graph)$name, network$nodes)]
  # Remove unwanted edge attributes
  if("seed" %in% edge_attr_names(network_graph)) network_graph <- delete_edge_attr(network_graph, "seed")
  if("ecc_rel_score" %in% edge_attr_names(network_graph)) network_graph <- delete_edge_attr(network_graph, "ecc_rel_score")

  if(is.null(node_table) == T) {
    node_table <- tibble("nodes" = as.character())
  }
  # similar to this to be done for edges when needed

  if(nrow(node_table) != 0) {
    for(i in 1:length(node_table)) {

      network_graph <- set_vertex_attr(graph = network_graph, name = colnames(node_table)[i],
                                       index = node_table$nodes, value = as.data.frame(node_table)[,i])
    }
  }

  network_filtered <- delete_vertices(network_graph, V(network_graph)$ecc_rel_score < relevance_score_cutoff) # needs to be universal if we have more scores
  }else {
    if(is.null(node_table) == T) {
      node_table <- tibble("nodes" = as.character())
    }
    # similar to this to be done for edges when needed

    if(nrow(node_table) != 0) {
      for(i in 1:length(node_table)) {

        network_graph <- set_vertex_attr(graph = network_graph, name = colnames(node_table)[i],
                                         index = node_table$nodes, value = as.data.frame(node_table)[,i])
      }
    }
    network_filtered <- network_graph
  }
  return(network_filtered)

}
