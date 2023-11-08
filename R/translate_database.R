#' Translate Database IDs
#'
#' This function translates and standardizes the given protein-protein interaction database ("inweb" or "string")
#' to ensure consistent usage of ID types across multiple databases. The databases are processed and enriched
#' with multiple identifiers like Ensembl Peptide ID, UniProt ID, Ensembl Gene ID, Ensembl Transcript ID, and HGNC symbols.
#'
#' @param database Either a character string specifying the type of database to translate - Accepted values are "inweb" or "string" - or a user-made
#' database with three columns, named "node", "interactor", "confidence_score".
#' @param id_type In case database is nor "string" or "inweb". A character stating which IDs are used in the submitted database.
#' Must be one among 'ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot'.
#'
#' @return A tibble with processed and translated databases enriched with multiple identifiers.
#'         The tibble contains columns with original identifiers as well as translated identifiers such as
#'         Ensembl Peptide IDs, UniProt IDs, Ensembl Gene IDs, Ensembl Transcript IDs, and HGNC symbols.
#' @export
#'
#' @importFrom dplyr mutate filter left_join rename as_tibble mutate_all na_if
#' @importFrom biomaRt useMart getBM
#' @importFrom stringr str_replace
#' @importFrom httr set_config config timeout
#' @import curl
#'
#' @examples
#' \dontrun{
#' # Example with the 'string' database
#' translated_string_db <- translate_database("string")
#'
#' # Example with the 'inweb' database
#' translated_inweb_db <- translate_database("inweb")
#'}
#' @details
#' The function performs several major tasks:
#' - Sets a timeout configuration for the fetching process.
#' - Validates the provided database type.
#' - Processes the chosen database by translating and standardizing identifiers.
#' - Handles missing and duplicated identifiers.
#' - Joins converted IDs to the original database.
#'
#' The function might encounter issues such as missing IDs or duplicated mappings.
#' IDs that couldn’t be mapped to the database and IDs that returned multiple mappings
#' will be reported to help in troubleshooting and data cleaning.
translate_database <- function(database, id_type) {

  # Check if database is a character string
  if (is.character(database)) {
    if (!database %in% c("inweb", "string")) {
      stop("The provided database name is not supported: choose either 'inweb' or 'string'.")
    }
    # If the database is a character and valid ('inweb' or 'string')
    if(database == "string") {
      # string <- read_table("../VirtualPulldownR/data/string_9606proteinlinksv12.0.txt")
      # remove the 9606 string from the entry-names and select the two column describing the interactions.
      # note: the "combined score" column contains confidence scores. The scores go from 1 to 1000 (they are mutiplied by 1000 to make them integers).
      # According to the STRING website (https://string-db.org/cgi/help?sessionId=bpCLNesEZjq4), proteins with < 150 have low confidence, medium confidence is 400,
      # high confidence is 700 and proteins with highest confidence score have > 900.

      string <- string %>%
        mutate(
          protein1 = str_replace(protein1, "9606\\.", ""),
          protein2 = str_replace(protein2, "9606\\.", "")
        )

      colnames(string)[3] <- "string_cs"

      ###
      # convert entries into multiple IDs
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

      timeout(seconds = 500)
      # h <- new_handle()
      # handle_setopt(h, timeout_ms = 500000)  # sets timeout to 500 seconds
      #options(timeout = 500)
      #set_config(config(timeout = 500))  # Set timeout to 500 seconds for http
      converted_ids_string1 <- getBM(
        attributes = c('ensembl_peptide_id', 'uniprotswissprot', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'),
        filters = 'ensembl_peptide_id',
        values = unique(c(string$protein1)),
        mart = ensembl
      )

      # I am losing some IDs because some ensp IDs might not have corresponding uniprot IDs. Moreover, some Ensembl Peptide IDs
      # from string$protein1 might not be present in the current version of the Ensembl database I am querying. It is important to notice
      # that multiple mapping are also possible. A sigle ENSP might correspond to multiple UniProt (or ensg etc) entries.

      # check for missing IDs
      queried_ids1 <- unique(c(string$protein1))
      retrieved_ids1 <- converted_ids_string1$ensembl_peptide_id
      missing_ids1 <- setdiff(queried_ids1, retrieved_ids1)

      # check for duplicate mappings
      duplicated_ids1 <- converted_ids_string1$ensembl_peptide_id[duplicated(converted_ids_string1$ensembl_peptide_id)]

      converted_ids_string2 <- getBM(
        attributes = c('ensembl_peptide_id', 'uniprotswissprot', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'),
        filters = 'ensembl_peptide_id',
        values = unique(c(string$protein2)),
        mart = ensembl
      )

      # check for missing IDs
      queried_ids2 <- unique(c(string$protein2))
      retrieved_ids2 <- converted_ids_string2$ensembl_peptide_id
      missing_ids2 <- setdiff(queried_ids2, retrieved_ids2)

      # check for duplicate mappings
      duplicated_ids2 <- converted_ids_string2$ensembl_peptide_id[duplicated(converted_ids_string2$ensembl_peptide_id)]
      tot_missing <- unique(c(missing_ids1, missing_ids2))
      tot_duplicate <- unique(c(duplicated_ids1, duplicated_ids2))

      print(paste(length(tot_missing), "IDs could not be mapped to the database", sep = " "))
      print(paste(length(tot_duplicate), "IDs returned multiple mappings", sep = " "))

      # converting empty cells in NAs
      converted_ids_string1 <- converted_ids_string1 %>%
        mutate_all(~na_if(., ""))

      converted_ids_string2 <- converted_ids_string2 %>%
        mutate_all(~na_if(., ""))

      ### join the database of the converted IDs to the original string database
      string_all_ids <- string %>%
        left_join(converted_ids_string1, by = c("protein1" = "ensembl_peptide_id"), relationship = "many-to-many") %>%
        rename(nodes_ensp = protein1) %>%
        rename(nodes_uniprot = uniprotswissprot) %>%
        rename(nodes_ensg = ensembl_gene_id) %>%
        rename(nodes_enst = ensembl_transcript_id) %>%
        rename(nodes_hgnc = hgnc_symbol) %>%
        left_join(converted_ids_string2, by = c("protein2" = "ensembl_peptide_id"), relationship = "many-to-many") %>%
        rename(interactors_ensp = protein2) %>%
        rename(interactors_uniprot = uniprotswissprot) %>%
        rename(interactors_ensg = ensembl_gene_id) %>%
        rename(interactors_enst = ensembl_transcript_id) %>%
        rename(interactors_hgnc = hgnc_symbol)

      database <- string_all_ids

    }
    else if(database == "inweb") {
      #data("inweb")
      colnames(inweb)[3] <- "inweb_cs"

      ###
      # convert entries into multiple IDs
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

      converted_ids_inweb1 <- getBM(
        attributes = c('uniprotswissprot', 'ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'),
        filters = 'uniprotswissprot',
        values = unique(c(inweb$node)),
        mart = ensembl
      )

      # I am losing some IDs because some ensp IDs might not have corresponding uniprot IDs. Moreover, some Ensembl Peptide IDs
      # from string$protein1 might not be present in the current version of the Ensembl database I am querying. It is important to notice
      # that multiple mapping are also possible. A sigle ENSP might correspond to multiple UniProt (or ensg etc) entries.

      # check for missing IDs
      queried_ids1 <- unique(c(inweb$node))
      retrieved_ids1 <- converted_ids_inweb1$uniprotswissprot
      missing_ids1 <- setdiff(queried_ids1, retrieved_ids1)

      # check for duplicate mappings
      duplicated_ids1 <- converted_ids_inweb1$uniprotswissprot[duplicated(converted_ids_inweb1$uniprotswissprot)]

      converted_ids_inweb2 <- getBM(
        attributes = c('uniprotswissprot', 'ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'),
        filters = 'uniprotswissprot',
        values = unique(c(inweb$interactor)),
        mart = ensembl
      )

      # check for missing IDs
      queried_ids2 <- unique(c(inweb$interactor))
      retrieved_ids2 <- converted_ids_inweb2$uniprotswissprot
      missing_ids2 <- setdiff(queried_ids2, retrieved_ids2)

      # check for duplicate mappings
      duplicated_ids2 <- converted_ids_inweb2$uniprotswissprot[duplicated(converted_ids_inweb2$uniprotswissprot)]
      tot_missing <- unique(c(missing_ids1, missing_ids2))
      tot_duplicate <- unique(c(duplicated_ids1, duplicated_ids2))

      print(paste(length(tot_missing), "IDs could not be mapped to the database", sep = " "))
      print(paste(length(tot_duplicate), "IDs returned multiple mappings", sep = " "))

      # converting empty cells in NAs
      converted_ids_inweb1 <- converted_ids_inweb1 %>%
        mutate_all(~na_if(., ""))

      converted_ids_inweb2 <- converted_ids_inweb2 %>%
        mutate_all(~na_if(., ""))

      ### join the database of the converted IDs to the original string database
      inweb_all_ids <- inweb %>%
        left_join(converted_ids_inweb1, by = c("node" = "uniprotswissprot"), relationship = "many-to-many") %>%
        rename(nodes_uniprot = node) %>%
        rename(nodes_ensp = ensembl_peptide_id) %>%
        rename(nodes_ensg = ensembl_gene_id) %>%
        rename(nodes_enst = ensembl_transcript_id) %>%
        rename(nodes_hgnc = hgnc_symbol) %>%
        left_join(converted_ids_inweb2, by = c("interactor" = "uniprotswissprot"), relationship = "many-to-many") %>%
        rename(interactors_uniprot = interactor) %>%
        rename(interactors_ensp = ensembl_peptide_id) %>%
        rename(interactors_ensg = ensembl_gene_id) %>%
        rename(interactors_enst = ensembl_transcript_id) %>%
        rename(interactors_hgnc = hgnc_symbol)

      database <- inweb_all_ids

    }

  } else if (is.data.frame(database)) {
    allowed_id_types <- c('ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot')

    # Check if id_type is a character string and one of the allowed values
    if (!is.character(id_type) || length(id_type) != 1 || !(id_type %in% allowed_id_types)) {
      stop(sprintf("id_type must be a single character string and one of the following: %s",
                   paste(shQuote(allowed_id_types), collapse = ", ")))
    }
    # check for required columns
    required_columns <- c("node", "interactor", "confidence_score")
    if (!all(required_columns %in% names(database))) {
      stop("The provided database does not have the required columns: 'node', 'interactor', 'confidence_score'.")
    }

    # convert entries into multiple IDs
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    timeout(seconds = 500)
    # h <- new_handle()
    # handle_setopt(h, timeout_ms = 500000)  # sets timeout to 500 seconds
    #options(timeout = 500)
    #set_config(config(timeout = 500))  # Set timeout to 500 seconds for http
    converted_ids1 <- getBM(
      attributes = c(id_type, 'ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot'),
      filters = id_type,
      values = unique(c(database$node)),
      mart = ensembl
    )
    # remove the duplicate column
    converted_ids1 <- converted_ids1 %>%
      select(-ends_with(".1"))

    # I am losing some IDs because some ensp IDs might not have corresponding uniprot IDs. Moreover, some Ensembl Peptide IDs
    # from string$protein1 might not be present in the current version of the Ensembl database I am querying. It is important to notice
    # that multiple mapping are also possible. A sigle ENSP might correspond to multiple UniProt (or ensg etc) entries.

    # check for missing IDs
    queried_ids1 <- unique(c(database$node))
    retrieved_ids1 <- converted_ids1[[id_type]]
    missing_ids1 <- setdiff(queried_ids1, retrieved_ids1)

    # check for duplicate mappings
    duplicated_ids1 <- converted_ids1[[id_type]][duplicated(converted_ids1[[id_type]])]

    timeout(seconds = 500)
    # h <- new_handle()
    # handle_setopt(h, timeout_ms = 500000)  # sets timeout to 500 seconds
    #options(timeout = 500)
    #set_config(config(timeout = 500))  # Set timeout to 500 seconds for http
    converted_ids2 <- getBM(
      attributes = c(id_type, 'ensembl_peptide_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot'),
      filters = id_type,
      values = unique(c(database$interactor)),
      mart = ensembl
    )
    # remove the duplicate column
    converted_ids2 <- converted_ids2 %>%
      select(-ends_with(".1"))

    # check for missing IDs
    queried_ids2 <- unique(c(database$interactor))
    retrieved_ids2 <- converted_ids2[[id_type]]
    missing_ids2 <- setdiff(queried_ids2, retrieved_ids2)

    # check for duplicate mappings
    duplicated_ids2 <- converted_ids2[[id_type]][duplicated(converted_ids2[[id_type]])]

    tot_missing <- unique(c(missing_ids1, missing_ids2))
    tot_duplicate <- unique(c(duplicated_ids1, duplicated_ids2))

    print(paste(length(tot_missing), "IDs could not be mapped to the database", sep = " "))
    print(paste(length(tot_duplicate), "IDs returned multiple mappings", sep = " "))

    # converting empty cells in NAs
    converted_ids1 <- converted_ids1 %>%
      mutate_all(~na_if(., ""))

    converted_ids2 <- converted_ids2 %>%
      mutate_all(~na_if(., ""))

    ### join the database of the converted IDs to the original string database
    all_ids <- database %>%
      left_join(converted_ids1, by = c("node" = id_type), relationship = "many-to-many") %>%
      rename(nodes_uniprot = node) %>%
      rename(nodes_ensp = ensembl_peptide_id) %>%
      rename(nodes_ensg = ensembl_gene_id) %>%
      rename(nodes_enst = ensembl_transcript_id) %>%
      rename(nodes_hgnc = hgnc_symbol) %>%
      left_join(converted_ids2, by = c("interactor" = id_type), relationship = "many-to-many") %>%
      rename(interactors_uniprot = interactor) %>%
      rename(interactors_ensp = ensembl_peptide_id) %>%
      rename(interactors_ensg = ensembl_gene_id) %>%
      rename(interactors_enst = ensembl_transcript_id) %>%
      rename(interactors_hgnc = hgnc_symbol)

    database <- all_ids

  } else {
    stop("The 'database' argument must be a character string specifying a supported database or a data frame with the required columns.")
  }

  # Return the processed database
  return(database)
}
