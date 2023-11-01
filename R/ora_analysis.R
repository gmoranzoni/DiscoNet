#' Convert Various Gene ID Types to Entrez IDs
#'
#' This function facilitates the conversion of various types of gene IDs (like "uniprot", "ensembl", etc.) to Entrez IDs.
#' It leverages the biomaRt package to achieve this.
#'
#' @param gene_list A character vector representing a list of gene identifiers to be converted.
#' @param mart A biomaRt object that specifies the database and dataset to use for the conversion.
#' @param id_type A character string that indicates the type of gene ID provided in `gene_list`
#'                (e.g., "uniprot", "ensembl", "ensembl_peptide", "ensembl_transcript", "hgnc").
#'
#' @return A character vector of converted Entrez gene IDs.
#' @importFrom biomaRt getBM
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' # Assuming `sample_genes` is a predefined list of gene identifiers and `sample_mart` is a defined biomaRt object:
#' # entrez_ids <- convert_to_entrez(gene_list = sample_genes, mart = sample_mart, id_type = "uniprot")
#' }
#'
convert_to_entrez <- function(gene_list, mart, id_type) {

  # Define a lookup list for the appropriate filters based on the ID type
  id_filters <- list(
    "uniprot" = "uniprot_gn_id",
    "ensembl" = "ensembl_gene_id",
    "ensembl_peptide" = "ensembl_peptide_id",
    "ensembl_transcript" = "ensembl_transcript_id",
    "hgnc" = "hgnc_symbol"
  )

  # Check if id_type is valid
  if(!id_type %in% names(id_filters)) {
    stop(paste("Invalid id_type:", id_type, ". Valid types are:", paste(names(id_filters), collapse=", ")))
  }

  res <- getBM(attributes = c("entrezgene_id", id_filters[[id_type]]),
               filters = id_filters[[id_type]],
               values = gene_list,
               mart = mart,
               uniqueRows = T)

  return(res$entrezgene_id %>% unique())
}


#' Over-Representation Analysis (ORA) on Network Communities
#'
#' Performs over-representation analysis (ORA) on communities derived from a network graph.
#' It supports various types of gene IDs and can fetch gene sets from the Molecular Signatures Database (MSigDB)
#' or use custom gene sets.
#'
#' @param graph An igraph object representing the network.
#' @param communities A list of igraph subgraphs representing the network communities.
#' @param organism A string indicating the organism, e.g., "human" or "mouse".
#' @param collection A string indicating the MSigDB collection, or NULL if not using MSigDB.
#' @param ontology A string indicating the MSigDB ontology, or NULL if not using this ontology.
#' @param custom_gene_set A logical flag, TRUE if using a custom gene set, FALSE otherwise. Defaults to FALSE.
#' @param universe A string indicating the universe for the ORA. Defaults to "whole_network".
#' @param id_type A string indicating the type of gene ID used, e.g., "uniprot", "ensembl", etc.
#' @import igraph
#' @importFrom biomaRt useMart
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr %>%
#' @importFrom purrr map
#' @importFrom fgsea fora
#' @return A list containing the ORA results for each community.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using example igraph object for graph and dummy communities.
#' # This is a placeholder; replace with real examples.
#'
#' # Assuming `sample_graph` is a predefined igraph object and `sample_communities` is a list of igraph subgraphs.
#' results <- ora_analysis(graph = sample_graph, communities = sample_communities, organism = "human", id_type = "uniprot")
#'}
ora_analysis <- function(graph, communities, organism, collection = NULL, ontology = NULL, custom_gene_set = F, universe = "whole_network", id_type) {
  stopifnot(is.igraph(graph))
  # Define a lookup list for the appropriate filters based on the ID type
  id_filters <- list(
    "uniprot" = "uniprot_gn_id",
    "ensembl" = "ensembl_gene_id",
    "ensembl_peptide" = "ensembl_peptide_id",
    "ensembl_transcript" = "ensembl_transcript_id",
    "hgnc" = "hgnc_symbol",
    "entrez" = "entrezgene_id"
  )

  # Check if id_type is valid
  if(!id_type %in% names(id_filters)) {
    stop(paste("Invalid id_type:", id_type, ". Valid types are:", paste(names(id_filters), collapse=", ")))
  }

  if(id_type != "entrez") {
    # Setting up biomaRt for ID conversion based on organism
    mart <- switch(organism,
                   "human" = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                   "mouse" = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl"),
                   stop(paste("Unsupported organism:", organism)))

    # Convert the universe to Entrez IDs
    if(universe == "whole_network") {
      universe_genes <- V(graph)$name
      universe <- convert_to_entrez(universe_genes, mart, id_type)
    }
  }

  # Check for the algorithm used
  # if length(communities) == 2, the algorithm used was mcode.
  if(length(communities) == 2) {
    communities <- communities$communities
  }

  # Convert community genes to Entrez IDs
  if(id_type != "entrez") {
    genes_subclusters <- lapply(communities, function(g) {
      genes <- as.character(V(g)$name)
      return(convert_to_entrez(genes, mart, id_type))
    })
  } else {
    genes_subclusters <- lapply(communities, function(g) {
      return(as.character(V(g)$name))
    })
  }

  if(custom_gene_set == F) {
    # Getting the gene_sets, from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
    gene_sets <- msigdbr(species = organism, category = collection, subcategory = ontology)
    # Convert the object into the correct format for fora
    gene_sets <- gene_sets %>% split(x = .$entrez_gene, f = .$gs_name)
    # Perform ORA
    results_ora <- genes_subclusters %>%
      map(., function(x) fora(genes = x, universe = universe, pathway = gene_sets))
  } else {
    # add code for custom gene sets if we want them
    results_ora <- NULL
  }

  return(results_ora)
}
