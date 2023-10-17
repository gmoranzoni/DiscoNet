#' Calculate ECC Relevance Scores for Network Nodes
#'
#' Calculates the ECC relevance score for all nodes in the provided network based on interactions
#' from a database and a set of seed nodes. This function considers triangle interactions
#' and combines interaction counts from the network and a provided database to determine
#' the ECC relevance scores.
#'
#' @param network A data frame representing the interactions in the network. Expected
#' to have columns 'nodes' and 'interactors' denoting the nodes and their interactions, respectively.
#' @param seed_nodes A vector containing seed nodes which are to be used in determining relevance scores.
#' @param database A data frame representing the interactions in the database. Expected
#' to have a column named 'nodes' representing the nodes.
#'
#' @return A data frame with nodes and their respective ECC relevance scores.
#'
#' @importFrom dplyr inner_join rename filter group_by summarise left_join rowwise mutate select
#' @importFrom tidyr replace_na
#'
#'
#' @examples
#' \dontrun{
#' # Assuming you have a data frame 'my_network' representing your network,
#' # a vector 'my_seed_nodes' containing seed nodes, and a data frame 'my_database'
#' # representing the interactions database:
#' scores <- calculate_relevance_score(my_network, my_seed_nodes, my_database)
#'}
calculate_relevance_score <- function(network, seed_nodes, database) {

  #nodes that share neighbor with a neighbor (in a triangle)
  triangle_interactions <- network %>%
    inner_join(., network, by="interactors", multiple = "all", relationship = "many-to-many") %>%
    rename(
      "nodes" = "nodes.x",
      "neighbors" = "interactors",
      "second_neighbors" = "nodes.y"
    ) %>%
    filter(nodes!=second_neighbors)  %>% #ensure that second neighbor is not self
    inner_join(., network, by=c("nodes"="nodes", "second_neighbors"="interactors"), multiple = "all", relationship = "many-to-many") %>%
    group_by(nodes, neighbors) %>%
    summarise(c = n(), .groups = "keep")

  # find total database number of interactions per node
  # I am considering the database of converted IDs
  database_interactions <- database %>%
    group_by(nodes) %>%
    summarise(int_total = n())

  # combine interaction count tables
  n_interactions <- network %>%
    left_join(., database_interactions, "nodes", multiple = "all")  %>%
    left_join(., database_interactions, by=c("interactors"="nodes"), multiple = "all")  %>%
    rowwise()  %>%
    mutate( int_total = min(int_total.x, int_total.y) ) %>%
    left_join(., triangle_interactions, by=c("nodes"="nodes", "interactors"="neighbors")) %>%
    replace_na(., replace = list("c" = 0)) %>%  #remove nan
    mutate(parts = c/int_total) %>%
    replace_na(., replace = list("parts" = 0)) %>%  #remove nan
    group_by(nodes) %>%
    summarise( score = mean(parts) )

  # calculate relevance score
  relevance_scores <- n_interactions %>% rowwise(.)  %>%
    mutate(seed = ifelse(nodes %in% seed_nodes == T, T, F)) %>% # test if seed nodes
    mutate("ecc_rel_score" = case_when(seed == T ~ 1,
                                       seed == F ~ score )) %>%
    select("nodes", "ecc_rel_score")

  return(relevance_scores)
}

#' Retrieve All Interactions from a Network Database
#'
#' This function retrieves all the interactions from the provided network database.
#' Since some networks might include one-way interactions, this function ensures
#' that all interactions are considered to be two-way (bidirectional) by adding
#' reverse interactions if they are not already present.
#'
#' @param database A data frame representing the interactions in the network.
#' Expected to have columns 'nodes' and 'interactors' denoting the nodes and their interactions, respectively.
#'
#' @return A data frame that includes all the interactions from the network, considering
#' both directions of each interaction. The returned data frame has columns 'nodes' and 'interactors'.
#'
#' @importFrom dplyr rename add_row distinct
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' # Assuming you have a data frame 'my_database' representing your network:
#' all_interactions <- get_all_interactions(my_database)
#'}
get_all_interactions <- function(database) {
  # The network might include one-way interactions.
  # When calculating relevance scores we assume all interactions to be two-way
  # Here we expand the network to include both directions of all interactions
  reverse_network <- database %>%
    rename(
      "nodes" = "interactors",
      "interactors" = "nodes"
    )

  network_all_interactions <- database %>%
    add_row(reverse_network) %>% # lengthen by reverse network
    distinct() %>% # remove duplicate rows
    na.omit()# remove NA

  return(network_all_interactions)
}

#' Conduct Virtual Pull-Down Analysis on a Protein Network
#'
#' This function conducts a virtual pull-down analysis on a protein-protein interaction network.
#' It fetches interactions from the given database (either 'inweb' or 'string'), converts them to the desired ID type,
#' and then identifies all possible interactions for the provided seed nodes. Additionally, it calculates relevance scores
#' for the identified interactions.
#'
#' @param seed_nodes A vector of seed proteins.
#' @param database A preprocessed database for analysis, either derived from the
#'        `translate_database` function or formatted similarly.
#' @param id_type The desired ID type for proteins. Can be one of: 'uniprot', 'ensg', 'ensp', 'enst', or 'hgnc'.
#' @param string_confidence_score A confidence score threshold for STRING database (default is 700). Should be set to NULL if not used.
#' @param zs_confidence_score ZS confidence score to use in the inweb database. Interactions with a score higher than 0.156 (default) are considered high confidence. Should be set to NULL if not used. Goes from 0 to 1.
#' @param order Numeric value (0 or 1) defining whether to fetch immediate interactors (order=1) or to also include second-level interactors (order=0).
#' @return A tibble that includes nodes, their interactors, calculated relevance scores, and a flag indicating
#' whether a node is a seed (1 for seed, 0 otherwise). Rows with NA values are removed from the final output,
#' and a message is displayed indicating the nodes that were removed.
#'
#' @importFrom dplyr %>% mutate filter select rename left_join full_join as_tibble mutate_all pull
#' @importFrom stats na.omit
#' @importFrom tidyr drop_na
#' @importFrom stringr str_replace
#' @importFrom biomaRt useMart getBM
#' @importFrom tibble tibble as_tibble as_tibble_col
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a vector of seed nodes:
#' seeds <- c("P12345", "Q67890", ...)
#' result <- virtual_pulldown(seed_nodes = seeds, database = "string", id_type = "uniprot")
#'}
virtual_pulldown <- function(seed_nodes, database, id_type, string_confidence_score = 700, zs_confidence_score = 0.156, order = 1) {

  if(!(order == 0 || order == 1)) {
    stop("Error: 'order' can only be either 0 or 1.")
  }

  id_filters <- c("uniprot", "ensg", "ensp", "enst", "hgnc")

  # Check if id_type is valid
  if(!id_type %in% id_filters) {
    stop(paste("Invalid id_type:", id_type, ". Valid types are:", paste(id_filters, collapse=", ")))
  }

  # filtering step based on the confidence score
  if("string_cs" %in% colnames(database)) {
    database <- database %>%
      filter(string_cs > string_confidence_score)
  }
  else if("inweb_cs" %in% colnames(database)) {
    database <- database %>%
      filter(inweb_cs > zs_confidence_score)
  }
  else {
    if(!is.null(zs_confidence_score) || !is.null(string_confidence_score)) {
      message("Neither 'string_cs' nor 'inweb_cs' cutoffs are provided. No filtering based on the confidence score is performed.")
    }
  }

  # select the relevant columns based on the IDs
  if(id_type == "ensp") {
    database <- database %>% select(nodes_ensp, interactors_ensp)
    database <- database %>% rename(
      "nodes" = "nodes_ensp",
      "interactors" = "interactors_ensp"
    )

  } else if(id_type == "uniprot") {
    database <- database %>% select(nodes_uniprot, interactors_uniprot)
    database <- database %>% rename(
      "nodes" = "nodes_uniprot",
      "interactors" = "interactors_uniprot"
    )

  } else if(id_type == "ensg") {
    database <- database %>% select(nodes_ensg, interactors_ensg)
    database <- database %>% rename(
      "nodes" = "nodes_ensg",
      "interactors" = "interactors_ensg"
    )

  } else if(id_type == "enst") {
    database <- database %>% select(nodes_enst, interactors_enst)
    database <- database %>% rename(
      "nodes" = "nodes_enst",
      "interactors" = "interactors_enst"
    )

  } else if(id_type == "hgnc") {
    database <- database %>% select(nodes_hgnc, interactors_hgnc)
    database <- database %>% rename(
      "nodes" = "nodes_hgnc",
      "interactors" = "interactors_hgnc"
    )
  }
  database <- database %>% na.omit

  if(order == 0) {

    network_final <- database %>%
      filter(nodes %in% seed_nodes & interactors %in% seed_nodes)

  }else if(order == 1) {
  all_network <- get_all_interactions(database = database)

  #### now the virtual pulldown itself

  seed_tibble <- tibble(nodes = seed_nodes, seed = 1)

  # find the interactors of the seed proteins in the database.
  # intermediate contains the seeds + the first order interactors of the seeds.
  intermediate <- all_network %>%
    filter((nodes %in% seed_tibble$nodes)|(interactors %in% seed_tibble$nodes)) %>%
    as_tibble()

  # all_proteins contains the whole list of proteins in the resulting network.
  all_proteins <- as_tibble_col(c(intermediate$nodes, intermediate$interactors),
                                column_name = "protein_name") %>% unique()

  # now, we want to find all the possible interactions in our network. So we need the 1st order interactors
  # of the seeds, and then the interactions happening between them.

  network <- all_network %>%
    filter(
      (nodes %in% all_proteins$protein_name) &
        (interactors %in% all_proteins$protein_name)
    ) %>%
    as_tibble()

  rel <- calculate_relevance_score(network = network, seed_nodes = seed_nodes, database = database)

  network_rel_score <- network %>%
    left_join(rel, by = "nodes") %>%
    full_join(seed_tibble, by = "nodes") %>%
    mutate(seed = replace_na(seed, 0)) # replace NAs with "0", only in the seed column. So that the non-seed are 0s.


  # Identify rows with NA and capture nodes
  na_nodes <- network_rel_score %>% filter(rowSums(is.na(.)) > 0) %>% pull(nodes)
  # Calculate number of rows with NA
  na_rows_count <- length(na_nodes)
  # Check for NAs and remove rows
  network_rel_score <- network_rel_score %>% drop_na()


  if(na_rows_count > 0) {
    nodes_str <- paste(na_nodes, collapse = ", ")
    cat(paste("Detected", na_rows_count, "nodes with no interactions:", nodes_str, ". The corresponding rows have been removed from the network.\n"))
  }

  network_final <- network_rel_score
  }

  return(network_final)
}
