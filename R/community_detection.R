#' Calculate Vertex Weights Using MCODE Vertex Weighting
#'
#' This function calculates the weights of vertices in a given graph based on the
#' k-core decomposition method used in MCODE (Molecular Complex Detection). The weights
#' reflect the relative importance or centrality of each vertex within the network,
#' aiding in community detection and other graph analysis tasks.
#'
#' @param graph An `igraph` object representing the input network. The function checks
#'              whether the input is a valid `igraph` object.
#'
#' @return A numeric vector containing the weights of each vertex in the input graph.
#'         The order of weights corresponds to the order of vertices in the graph.
#'
#' @importFrom igraph induced_subgraph graph.coreness ecount vcount which_loop
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Creating a sample graph
#' g <- make_ring(10)
#'
#' # Calculating vertex weights
#' weights <- mcode.vertex.weighting(g)
#'
#' print(weights)
#' }
mcode.vertex.weighting<-function(graph){
  stopifnot(is.igraph(graph))

  # get neighborhood of a vertex
  get_neighbors <- function(graph, v) {
    neighbors(graph, v)
  }

  # Define weight calculation for a given vertex
  calculate_weight <- function(graph, vertex) {
    # Get neighbors of the vertex
    neigh <- get_neighbors(graph, vertex)

    # If there are no neighbors, return 0
    if(length(neigh) == 0) return(0)

    # Create subgraph using these neighbors
    subg <- induced_subgraph(graph, vids = neigh)

    # Calculate coreness of the subgraph
    core <- graph.coreness(subg)

    # If all coreness values are NA, return 0
    if(all(is.na(core))) return(0)

    k <- max(core, na.rm = T)

    # Get k-core subgraph
    kcore <- induced_subgraph(subg, vids = which(core == k))

    # the following is the calculation of the density.
    # if there are self-loops in the graph density calculated in one way, otherwise other way. In string and inweb database there are no self loops.
    # Here, the density is calculated as the number of actual edges divided by the total possible edges.The density
    # is then multiplied by the k-core number k, which results in the weight value of the vertex.

    if (vcount(kcore) > 1) {
     # if (any(which_loop(kcore))) {
      #   return(k * ecount(kcore) / choose(vcount(kcore) + 1, 2))
      # } else {
        return(k * ecount(kcore) / choose(vcount(kcore), 2))
     # }
    } else {
      return(0)
    }
  }

  # Calculate weights for all vertices
  weights <- sapply(seq_len(vcount(graph)), function(vertex) {
    calculate_weight(graph, vertex)
  })

  return(weights)
}

#' MCODE Find Complex Procedure
#'
#' This function identifies a single densely connected cluster (or complex) within a graph
#' using the MCODE (Molecular Complex Detection) algorithm. A low-level helper function
#' that leverages a compiled C function to efficiently find regions of high connectivity.
#'
#' @param neigh A vector containing the neighboring vertices for each vertex in the graph.
#' @param neighbor_indx An index vector where each element points to the positions in
#'        the `neigh` vector, defining the range of neighbors for each vertex.
#' @param vertex_weight A numeric vector containing the weights for each vertex in the graph.
#' @param D A numeric value representing the density parameter utilized in the complex detection process.
#' @param seed_vertex An integer specifying the vertex where the search for the complex initiates.
#' @param seen A logical vector indicating whether each vertex has been processed or encountered previously.
#'
#' @useDynLib DiscoNet, .registration = TRUE
#' @return A list with two components:
#'         - `seen`: A logical vector updated to indicate which vertices have been processed.
#'         - `COMPLEX`: A numeric vector containing the indices of vertices included in the identified complex.
#'
#' @examples
#' \dontrun{
#' # Assuming pre-processed graph data: neigh, neighbor_indx, vertex_weight, D, seed_vertex, and seen are available.
#'
#' result <- mcode.find.complex(neigh, neighbor_indx, vertex_weight, D, seed_vertex, seen)
#' print(result$COMPLEX) # Output the vertices included in the identified complex
#' }
mcode.find.complex <- function(neigh, neighbor_indx, vertex_weight,
                               D, seed_vertex, seen) {

  # Initialize COMPLEX as a zero vector with the same length as seen
  COMPLEX <- integer(length(seen))

  res <- .C("mcode_complex",
            neighbor = as.integer(neigh),
            neighbor_indx = as.integer(neighbor_indx),
            vertex_weight = as.single(vertex_weight),
            D = as.single(D),
            seed_vertex = as.integer(seed_vertex),
            seen = as.integer(seen),
            COMPLEX = COMPLEX
  )

  return(list(seen=res$seen,COMPLEX=which(res$COMPLEX!=0)))
}

#' MCODE Find Complexes Procedure
#'
#' This function discovers multiple densely connected clusters (complexes) within
#' a given graph. It iterates over all vertices, utilizing the `mcode.find.complex`
#' helper function, to efficiently uncover all potential complexes based on vertex
#' weights and a specified density parameter.
#'
#' @param graph An `igraph` object representing the input network.
#' @param neigh A list where each element contains the neighbors of each vertex in the graph.
#' @param vertex_weight A numeric vector containing the weights for each vertex in the graph.
#' @param D A numeric value representing the density parameter used in the identification of complexes.
#'
#' @return A list with two components:
#'         - `COMPLEX`: A list where each element is a numeric vector of vertex indices forming an identified complex.
#'         - `seen`: A numeric vector marking the vertices that have been processed.
#'
#' @examples
#' \dontrun{
#' # Assuming that you have a graph and its related information like neigh, vertex_weight, and D.
#'
#' result <- mcode.find.complexes(graph, neigh, vertex_weight, D)
#' print(length(result$COMPLEX)) # Output the number of identified complexes
#' }
mcode.find.complexes<-function(graph, neigh, vertex_weight, D) {
  seen<-rep(0,vcount(graph))

  neigh<-lapply(neigh,function(item){item[-1]})
  neighbors.indx<-cumsum(unlist(lapply(neigh,length)))

  neighbors.indx<-c(0,neighbors.indx)
  neigh<-unlist(neigh)-1

  COMPLEX <- list()
  n <- 1
  w.order <- order(vertex_weight, decreasing = TRUE)

  for(i in w.order){
    if(!(seen[i])){
      res <- mcode.find.complex(neigh = neigh, neighbor_indx = neighbors.indx, vertex_weight, D, i-1, seen)
      if(length(res$COMPLEX)>1){
        COMPLEX[[n]]<-res$COMPLEX
        seen<-res$seen
        n<-n+1
      }
    }
  }
  rm(neigh)
  return(list(COMPLEX = COMPLEX,seen = seen))
}

#' MCODE Fluff Complex
#'
#' This function refines the identified complex resulting from the MCODE algorithm by considering
#' the neighborhood and density criteria. For each vertex in the identified complex, the function
#' examines its neighborhood. If the density of the neighborhood surpasses a specific threshold
#' (`fdt`), the neighboring vertices are incorporated into the complex.
#'
#' @param graph An `igraph` object representing the input graph.
#' @param vertex.weight A numeric vector containing weights for each vertex in the graph.
#' @param fdt A numeric threshold (default is 0.8) for neighborhood density. Neighbors with a density
#'           exceeding this threshold will be incorporated into the complex.
#' @param complex.g A numeric vector containing the vertex indices that form the identified complex.
#' @param seen A numeric vector indicating the vertices that have been processed or visited.
#'
#' @return A numeric vector containing the vertex indices that constitute the expanded complex.
#'
#' @importFrom igraph neighborhood induced.subgraph graph.density
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have an igraph object (graph), vertex weights (vertex.weight),
#' # an initial complex (complex.g), and seen vertices (seen).
#' fluffed_complex <- mcode.fluff.complex(graph, vertex.weight, fdt=0.8, complex.g, seen)
#' }
mcode.fluff.complex<-function(graph,vertex.weight,fdt = 0.8,complex.g,seen){
  seq_complex.g <- seq_along(complex.g) # Preparing to iterate through all vertices in complex C
  for(i in seq_complex.g){ # for all u in C do
    node.neighbor <- unlist(neighborhood(graph,1,complex.g[i])) # Finding the neighbors of each vertex in the complex
    if(length(node.neighbor) > 1){ # Ensuring there are neighbors to consider
      subg <- induced.subgraph(graph,node.neighbor) # Creating subgraph with neighbors
      if(graph.density(subg, loops=FALSE) > fdt){ # if weight of u > d (interpreted as density > fluff density threshold)
        complex.g <- c(complex.g,node.neighbor) # add u to complex C
      }
    }
  }

  return(unique(complex.g))
}

#' MCODE Post-Processing Procedure
#'
#' This function refines the detected clusters from the MCODE algorithm. It filters out smaller
#' clusters (size <= 2) and, based on given parameters, may expand ('fluff') clusters or trim
#' singly-connected nodes ('haircut') from them. This post-processing aims to improve the
#' relevance and robustness of the detected clusters.
#'
#' @param graph An `igraph` object representing the input network.
#' @param vertex.weight A numeric vector containing weights assigned to each vertex in the graph.
#' @param haircut Logical; determines whether to trim singly-connected nodes from the clusters.
#' @param fluff Logical; determines whether to expand clusters by including bordering nodes.
#' @param fdt A numeric value setting the threshold for the fluff procedure; default is 0.8.
#' @param set.complex.g A list where each element is a numeric vector representing the indices of
#'                      vertices in a detected cluster.
#' @param seen A numeric vector indicating which vertices have been processed or visited.
#'
#' @return A list where each element is a numeric vector of vertex indices, representing refined
#'         clusters after post-processing steps have been applied.
#'
#' @importFrom igraph neighborhood induced.subgraph graph.density graph.coreness
#' @examples
#' \dontrun{
#' # Assuming the existence of a graph, vertex weights, and other necessary inputs:
#' # graph, vertex.weight, haircut, fluff, set.complex.g, and seen.
#' refined_clusters <- mcode.post.process(graph, vertex.weight, haircut, fluff, fdt=0.8, set.complex.g, seen)
#' print(length(refined_clusters)) # Output the number of refined clusters.
#' }
mcode.post.process<-function(graph,vertex.weight,haircut,fluff,fdt=0.8, set.complex.g,seen) {
  # Pseudocode: if c not 2-core then filter
  indx<-unlist(lapply(set.complex.g,
                      function(complex.g){
                        if(length(complex.g)<=2)
                          0
                        else
                          1
                      }
  ))
  set.complex.g<-set.complex.g[indx!=0]
  set.complex.g<-lapply(set.complex.g,
                        function(complex.g){
                          coreness<-graph.coreness(induced.subgraph(graph,complex.g))
                          if(fluff){
                            complex.g<-mcode.fluff.complex(graph,vertex.weight,fdt,complex.g,seen)
                            if(haircut){
                              ## coreness needs to be recalculated
                              coreness<-graph.coreness(induced.subgraph(graph,complex.g))
                              complex.g<-complex.g[coreness>1]
                            }
                          }else if(haircut){
                            complex.g<-complex.g[coreness>1]
                          }
                          return(complex.g)
                        })
  set.complex.g<-set.complex.g[lapply(set.complex.g,length)>2]
  return(set.complex.g)
}




#' Clusters a network using the MCODE (Molecular Complex Detection) algorithm.
#' MCODE identifies densely connected regions in large protein interaction networks
#' that may represent molecular complexes.
#'
#' @title MCODE Network Clustering
#' @description This function applies the MCODE clustering algorithm to identify
#'   molecular complexes in protein interaction networks. It includes options for vertex
#'   weighting, cluster expansion, and more.
#'
#' @param graph An igraph object representing the network to be clustered.
#' @param D Numeric. Defines the vertex weight percentage for scoring vertices
#'   during the clustering process. Defaults to 0.5.
#' @param haircut Logical. If TRUE, singly-connected nodes are removed from clusters. Defaults to FALSE.
#' @param fluff Logical. If TRUE, expands the clusters by one neighbor shell outwards. Defaults to FALSE.
#' @param fdt Numeric. Sets the cluster density cutoff for fluffing procedure. Defaults to 0.8.
#' @param loops Logical. Determines whether to include self-loops in the graph. Defaults to TRUE.
#'
#' @references Bader GD, Hogue CW. An automated method for finding molecular complexes
#'   in large protein interaction networks. BMC Bioinformatics. 2003 Jan 13;4(1):2.
#' @return A list containing clusters sorted by a computed score based on size and edge count.
#'   Each element of the list represents a cluster and contains vertex IDs belonging to that cluster.
#' @importFrom igraph is.igraph simplify neighborhood induced.subgraph is.loop ecount vcount
#' @seealso \code{\link{cluster}}
#' @examples
#' \dontrun{
#'   nlocal <- data.frame(c("DVL1", "DVL2", "DVL3"))
#'   net <- construction(input = nlocal, db = "HPRD", species = "human",
#'                       ID.type = "Gene symbol", hierarchy = 1)
#'   mcode_result <- mcode(net, D = 0.9, haircut = TRUE, fluff = TRUE, fdt = 0.1)
#' }

# This is the main function that orchestrates the entire MCODE clustering process.
# First, it checks the validity of the input graph and parameters.
# Computes vertex weights using the vertex weighting function.
# Finds clusters using the find.complexes function.
# Post-processes these clusters with the post.process function.
# Finally, it computes a score for each cluster (based on size and edge count) and returns the clusters sorted by this score.
mcode <- function(graph, D = 0.5, haircut = FALSE, fluff = FALSE, fdt = 0.8, loops = TRUE) {

  stopifnot(is.igraph(graph))
  if(D > 1 | D < 0){
    stop("D must be between 0 and 1")
  }
  if(!loops){
    graph <- simplify(graph, remove.multiple = FALSE, remove.loops = TRUE)
  }

  neigh <- neighborhood(graph, 1)

  W <- mcode.vertex.weighting(graph)
  res <- mcode.find.complexes(graph, neigh = neigh, vertex_weight = W, D = D)
  COMPLEX <- mcode.post.process(graph,vertex.weight = W,haircut = haircut, fluff = fluff,
                                fdt = fdt, res$COMPLEX, res$seen)
  score <- unlist(lapply(COMPLEX,
                         function(complex.g){
                           complex.g<-induced.subgraph(graph,complex.g)
                           if(any(is.loop(complex.g)))
                             score <- ecount(complex.g)/choose(vcount(complex.g) + 1,2) * vcount(complex.g)
                           else
                             score <- ecount(complex.g)/choose(vcount(complex.g), 2) * vcount(complex.g)
                           return(score)
                         }
  ))
  order_score <- order(score, decreasing = TRUE)
  return(list(COMPLEX = COMPLEX[order_score], score = score[order_score]))
}


# community detection, read manual, you can also consider weights == edge attributes.
# based on paper Topological and functional comparison of community detection algorithms in
# biological networks (Sara Rahiminejad et al.) we selected possible community detection
# algorithms that can be selected in our tool.

#' Community Detection in Graphs
#'
#' Performs community detection on a given graph using the specified algorithm.
#' Supports the Louvain, Leiden, and MCODE algorithms. For the MCODE algorithm,
#' additional post-processing options like `fluff` and `haircut` are available.
#'
#' @title Community Detection in Graphs
#' @description Applies selected community detection algorithm on a network and
#'   returns detected communities as induced subgraphs. Options are available for
#'   various algorithms and their respective parameters.
#'
#' @param network An `igraph` object. The input network for community detection.
#' @param algorithm Character string. Specifies the algorithm ('louvain', 'leiden', or 'mcode').
#' @param D Numeric. Vertex weight percentage, relevant for MCODE. Default is 0.5.
#' @param haircut Logical. MCODE-specific, trims singly-connected nodes. Default is FALSE.
#' @param fluff Logical. MCODE-specific, expands clusters to include bordering nodes. Default is FALSE.
#' @param fdt Numeric. MCODE-specific, used in the fluff procedure. Default is 0.8.
#' @param loops Logical. Indicates if loops are allowed in the graph. Default is FALSE.
#' @param resolution Numeric. Resolution parameter for algorithms like Louvain and Leiden,
#'   affecting the granularity of communities. Default is 1.
#'
#' @return A list of induced subgraphs representing detected communities. For MCODE,
#'   the list also includes scores for each subgraph.
#' @importFrom igraph is.igraph V cluster_louvain membership cluster_leiden induced_subgraph
#' @export
#' @examples
#' \dontrun{
#' # Assuming an igraph object 'my_graph':
#' result_louvain <- community_detection(my_graph, algorithm = "louvain")
#' result_leiden <- community_detection(my_graph, algorithm = "leiden")
#' result_mcode <- community_detection(my_graph, algorithm = "mcode", D = 0.6, fluff = TRUE)
#'}

community_detection <- function(network, algorithm, D = 0.5, haircut = FALSE, fluff = FALSE, fdt = 0.8, loops = FALSE, resolution = 1) {

  stopifnot(is.igraph(network))

  # Check if the provided algorithm is valid
  if(!algorithm %in% c("louvain", "leiden", "mcode")) {
    stop("The provided algorithm is not supported. Choose either 'mcode', 'louvain', or 'leiden'.")
  }

  if(algorithm == "louvain") {
    comm_det <- cluster_louvain(network, resolution = resolution)
    # Extracting membership information from communities
    mem <- membership(comm_det)
    # Generate a list of vertices for each community
    communities <- split(V(network), mem)

  } else if(algorithm == "leiden") {
    comm_det <- cluster_leiden(network, resolution_parameter = resolution)
    # Extracting membership information from communities
    mem <- membership(comm_det)
    # Generate a list of vertices for each community
    communities <- split(V(network), mem)

  } else if(algorithm == "mcode") {
    comm_det <- mcode(network, D = D, haircut = haircut, fluff = fluff, fdt = fdt, loops = loops)
    # Extracting the communities
    communities <- comm_det$COMPLEX
    scores <- comm_det$score
  } else {
    stop("Unknown algorithm selected")
  }

  if(algorithm != "mcode") {
    # Create induced subgraphs
    list_of_subgraph <- lapply(communities, function(community) {
      induced_subgraph(network, vids = community)
    })
  } else if(algorithm == "mcode") {
    list_of_subgraph <- list()
    # Create induced subgraphs
    list_of_subgraph$communities <- lapply(communities, function(community) {
      induced_subgraph(network, vids = community)
    })
    list_of_subgraph$scores <- scores
  }


  return(list_of_subgraph)

}
