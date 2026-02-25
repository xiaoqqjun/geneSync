#' Gene Symbol and ID Conversion
#'
#' Convert gene symbols, Ensembl IDs, or Entrez IDs to other formats with synonym support.
#' When converting symbols, the function first attempts exact matching
#' against official symbols, then falls back to synonym matching.
#'
#' @param genes Character vector of gene symbols, Ensembl IDs, or Entrez IDs to convert.
#' @param species Species identifier. Must be one of "homo" (human),
#'   "mmu" (mouse), or "rat". Default is "homo".
#' @param from Input format. Must be one of "symbol", "ensembl_id", or "gene_id".
#'   Default is "symbol".
#' @param to Output format. Must be one of "symbol", "ensembl_id", or "gene_id".
#'   Default is "symbol" (normalizes gene symbols).
#' @param direction (Deprecated) Conversion direction. "symbol_to_id" converts gene symbols
#'   to Entrez IDs, "id_to_symbol" converts Entrez IDs to official symbols.
#'   Use \code{from} and \code{to} parameters instead.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{input}{Original input gene name or ID}
#'   \item{symbol}{NCBI gene symbol}
#'   \item{symbol_authority}{Official symbol from nomenclature authority (HGNC/MGI/RGD), NA if not available}
#'   \item{gene_id}{Entrez Gene ID}
#'   \item{ensembl_id}{Ensembl gene ID (if available)}
#'   \item{chromosome}{Chromosome location}
#'   \item{match_type}{Type of match: "exact_authority", "exact_symbol", "synonym", "id_match", "ensembl_match", "not_found", or "invalid_id"}
#'   \item{final_symbol}{Recommended symbol: symbol_authority > symbol > input (fallback)}
#' }
#'
#' @details
#' The matching priority for symbol input:
#' \enumerate{
#'   \item Exact match against symbol_authority (case-insensitive)
#'   \item Exact match against NCBI symbol (case-insensitive)
#'   \item Match against Synonyms (case-insensitive)
#' }
#'
#' The \code{final_symbol} column provides the recommended symbol to use:
#' \itemize{
#'   \item If symbol_authority is available, use symbol_authority
#'   \item Otherwise, if symbol is available, use symbol
#'   \item Otherwise, keep the original input (for unmatched genes)
#' }
#'
#' @examples
#' # Convert human gene symbols to IDs
#' gene_convert(c("GAPDH", "TP53", "BRCA1"), species = "homo")
#'
#' # Mitochondrial gene: ATP6 -> MT-ATP6
#' gene_convert("ATP6", species = "homo")
#'
#' # Convert using synonyms
#' gene_convert(c("G3PD", "p53"), species = "homo")
#'
#' # Convert IDs to symbols (new from/to syntax)
#' gene_convert(c("2597", "7157"), species = "homo", from = "gene_id", to = "symbol")
#'
#' # Convert Ensembl IDs to symbols
#' gene_convert(c("ENSG00000111640", "ENSG00000141510"), species = "homo", from = "ensembl_id")
#'
#' # Convert symbols to Ensembl IDs
#' gene_convert(c("GAPDH", "TP53"), species = "homo", to = "ensembl_id")
#'
#' # Mouse genes
#' gene_convert(c("Gapdh", "Tp53"), species = "mmu")
#'
#' # Backward compatibility with direction parameter (deprecated)
#' gene_convert(c("GAPDH", "TP53"), species = "homo", direction = "symbol_to_id")
#'
#' @export
gene_convert <- function(genes,
                         species = c("homo", "mmu", "rat"),
                         from = c("symbol", "ensembl_id", "gene_id"),
                         to = c("symbol", "ensembl_id", "gene_id"),
                         direction = NULL) {

  # Parameter validation
  species <- match.arg(species)

  # Backward compatibility: handle direction parameter
  if (!is.null(direction)) {
    direction <- match.arg(direction, c("symbol_to_id", "id_to_symbol"))
    if (missing(from) && missing(to)) {
      if (direction == "symbol_to_id") {
        from <- "symbol"
        to <- "gene_id"
      } else {
        from <- "gene_id"
        to <- "symbol"
      }
    }
  }

  # Set defaults
  if (missing(from) || is.null(from)) {
    from <- "symbol"
  } else {
    from <- match.arg(from)
  }

  if (missing(to) || is.null(to)) {
    to <- "symbol"
  } else {
    to <- match.arg(to)
  }

  if (length(genes) == 0) {
    stop("Input 'genes' cannot be empty.", call. = FALSE)
  }

  # Load species-specific database
  gene_db <- switch(species,
    "homo" = gene_db_homo,
    "mmu"  = gene_db_mmu,
    "rat"  = gene_db_rat
  )

  # Convert to character and get unique values
  genes <- as.character(genes)
  unique_genes <- unique(genes)

  # Perform conversion based on input type
  if (from == "symbol") {
    result <- .symbol_to_id_fast(unique_genes, gene_db)
  } else if (from == "gene_id") {
    result <- .id_to_symbol_fast(unique_genes, gene_db)
  } else if (from == "ensembl_id") {
    result <- .ensembl_to_info_fast(unique_genes, gene_db)
  }

  # Add final_symbol column: symbol_authority > symbol > input
  result$final_symbol <- ifelse(
    !is.na(result$symbol_authority),
    result$symbol_authority,
    ifelse(!is.na(result$symbol),
           result$symbol,
           result$input)
  )

  # Extract requested output column if different from default
  # For to = "symbol", we already have final_symbol
  # For to = "gene_id" or "ensembl_id", the columns are already in result

  # Reset row names
  rownames(result) <- NULL

  return(result)
}


#' Fast Symbol to ID Conversion (Internal)
#' @param genes Character vector of gene symbols
#' @param gene_db Gene database data frame
#' @return Data frame with conversion results
#' @keywords internal
.symbol_to_id_fast <- function(genes, gene_db) {

  n_genes <- length(genes)

  # Pre-compute uppercase versions for matching (do once)
  genes_upper <- toupper(genes)
  db_symbol_upper <- toupper(gene_db$symbol)
  db_authority_upper <- toupper(gene_db$symbol_authority)
  db_synonyms_upper <- toupper(gene_db$synonyms)

  # Initialize result vectors
  result_symbol <- rep(NA_character_, n_genes)
  result_authority <- rep(NA_character_, n_genes)
  result_gene_id <- rep(NA_integer_, n_genes)
  result_ensembl <- rep(NA_character_, n_genes)
  result_chr <- rep(NA_character_, n_genes)
  result_match_type <- rep("not_found", n_genes)

  # Track which genes still need matching
  unmatched <- rep(TRUE, n_genes)

  # Step 1: Exact match against symbol_authority (vectorized using match)
  authority_lookup <- setNames(seq_len(nrow(gene_db)), db_authority_upper)
  authority_match <- authority_lookup[genes_upper]

  matched_auth <- !is.na(authority_match)
  if (any(matched_auth)) {
    idx <- authority_match[matched_auth]
    result_symbol[matched_auth] <- gene_db$symbol[idx]
    result_authority[matched_auth] <- gene_db$symbol_authority[idx]
    result_gene_id[matched_auth] <- gene_db$gene_id[idx]
    result_ensembl[matched_auth] <- gene_db$ensembl_id[idx]
    result_chr[matched_auth] <- gene_db$chromosome[idx]
    result_match_type[matched_auth] <- "exact_authority"
    unmatched[matched_auth] <- FALSE
  }

  # Step 2: Exact match against symbol (NCBI)
  if (any(unmatched)) {
    symbol_lookup <- setNames(seq_len(nrow(gene_db)), db_symbol_upper)
    symbol_match <- symbol_lookup[genes_upper]

    matched_sym <- !is.na(symbol_match) & unmatched
    if (any(matched_sym)) {
      idx <- symbol_match[matched_sym]
      result_symbol[matched_sym] <- gene_db$symbol[idx]
      result_authority[matched_sym] <- gene_db$symbol_authority[idx]
      result_gene_id[matched_sym] <- gene_db$gene_id[idx]
      result_ensembl[matched_sym] <- gene_db$ensembl_id[idx]
      result_chr[matched_sym] <- gene_db$chromosome[idx]
      result_match_type[matched_sym] <- "exact_symbol"
      unmatched[matched_sym] <- FALSE
    }
  }

  # Step 3: Match against Synonyms (only for remaining unmatched - slower, use loop)
  if (any(unmatched)) {
    for (i in which(unmatched)) {
      escaped_gene <- .escape_regex_fixed(genes_upper[i])
      pattern <- paste0("(^|\\|)", escaped_gene, "(\\||$)")
      match_idx <- grep(pattern, db_synonyms_upper)
      if (length(match_idx) > 0) {
        idx <- match_idx[1]
        result_symbol[i] <- gene_db$symbol[idx]
        result_authority[i] <- gene_db$symbol_authority[idx]
        result_gene_id[i] <- gene_db$gene_id[idx]
        result_ensembl[i] <- gene_db$ensembl_id[idx]
        result_chr[i] <- gene_db$chromosome[idx]
        result_match_type[i] <- "synonym"
        unmatched[i] <- FALSE
      }
    }
  }

  # Build result data frame
  data.frame(
    input = genes,
    symbol = result_symbol,
    symbol_authority = result_authority,
    gene_id = result_gene_id,
    ensembl_id = result_ensembl,
    chromosome = result_chr,
    match_type = result_match_type,
    stringsAsFactors = FALSE
  )
}


#' Fast ID to Symbol Conversion (Internal)
#' @param genes Character vector of gene IDs
#' @param gene_db Gene database data frame
#' @return Data frame with conversion results
#' @keywords internal
.id_to_symbol_fast <- function(genes, gene_db) {

  n_genes <- length(genes)

  # Convert to integers
  gene_ids <- suppressWarnings(as.integer(genes))

  # Initialize result vectors
  result_symbol <- rep(NA_character_, n_genes)
  result_authority <- rep(NA_character_, n_genes)
  result_gene_id <- rep(NA_integer_, n_genes)
  result_ensembl <- rep(NA_character_, n_genes)
  result_chr <- rep(NA_character_, n_genes)
  result_match_type <- rep("not_found", n_genes)

  # Mark invalid IDs
  invalid_idx <- is.na(gene_ids)
  result_match_type[invalid_idx] <- "invalid_id"

  # Match valid IDs using lookup
  valid_idx <- !invalid_idx
  if (any(valid_idx)) {
    id_lookup <- setNames(seq_len(nrow(gene_db)), gene_db$gene_id)
    id_match <- id_lookup[as.character(gene_ids[valid_idx])]

    matched <- !is.na(id_match)
    matched_positions <- which(valid_idx)[matched]

    if (length(matched_positions) > 0) {
      idx <- id_match[matched]
      result_symbol[matched_positions] <- gene_db$symbol[idx]
      result_authority[matched_positions] <- gene_db$symbol_authority[idx]
      result_gene_id[matched_positions] <- gene_db$gene_id[idx]
      result_ensembl[matched_positions] <- gene_db$ensembl_id[idx]
      result_chr[matched_positions] <- gene_db$chromosome[idx]
      result_match_type[matched_positions] <- "id_match"
    }
  }

  # Build result data frame
  data.frame(
    input = genes,
    symbol = result_symbol,
    symbol_authority = result_authority,
    gene_id = result_gene_id,
    ensembl_id = result_ensembl,
    chromosome = result_chr,
    match_type = result_match_type,
    stringsAsFactors = FALSE
  )
}


#' Fast Ensembl ID to Info Conversion (Internal)
#' @param genes Character vector of Ensembl IDs
#' @param gene_db Gene database data frame
#' @return Data frame with conversion results
#' @keywords internal
.ensembl_to_info_fast <- function(genes, gene_db) {

  n_genes <- length(genes)

  # Normalize Ensembl IDs: remove version suffix and convert to uppercase
  # ENSG00000139618.9 -> ENSG00000139618
  genes_normalized <- toupper(genes)
  genes_normalized <- sub("\\.\\d+$", "", genes_normalized)  # Remove version suffix

  # Also normalize database Ensembl IDs
  db_ensembl_normalized <- toupper(gene_db$ensembl_id)
  db_ensembl_normalized <- sub("\\.\\d+$", "", db_ensembl_normalized)

  # Initialize result vectors
  result_symbol <- rep(NA_character_, n_genes)
  result_authority <- rep(NA_character_, n_genes)
  result_gene_id <- rep(NA_integer_, n_genes)
  result_ensembl <- rep(NA_character_, n_genes)
  result_chr <- rep(NA_character_, n_genes)
  result_match_type <- rep("not_found", n_genes)

  # Mark invalid Ensembl IDs (basic pattern check)
  # Valid patterns: ENSG..., ENSMUSG..., ENSRNOG...
  valid_ensembl_pattern <- "^ENS[A-Z]*G\\d+$"
  invalid_idx <- !grepl(valid_ensembl_pattern, genes_normalized)
  result_match_type[invalid_idx] <- "invalid_ensembl"

  # Match valid Ensembl IDs using lookup
  valid_idx <- !invalid_idx
  if (any(valid_idx)) {
    # Create lookup table (skip NA values in database)
    valid_db_idx <- !is.na(db_ensembl_normalized) & db_ensembl_normalized != ""
    ensembl_lookup <- setNames(
      which(valid_db_idx),
      db_ensembl_normalized[valid_db_idx]
    )

    ensembl_match <- ensembl_lookup[genes_normalized[valid_idx]]

    matched <- !is.na(ensembl_match)
    matched_positions <- which(valid_idx)[matched]

    if (length(matched_positions) > 0) {
      idx <- ensembl_match[matched]
      result_symbol[matched_positions] <- gene_db$symbol[idx]
      result_authority[matched_positions] <- gene_db$symbol_authority[idx]
      result_gene_id[matched_positions] <- gene_db$gene_id[idx]
      result_ensembl[matched_positions] <- gene_db$ensembl_id[idx]
      result_chr[matched_positions] <- gene_db$chromosome[idx]
      result_match_type[matched_positions] <- "ensembl_match"
    }
  }

  # Build result data frame
  data.frame(
    input = genes,
    symbol = result_symbol,
    symbol_authority = result_authority,
    gene_id = result_gene_id,
    ensembl_id = result_ensembl,
    chromosome = result_chr,
    match_type = result_match_type,
    stringsAsFactors = FALSE
  )
}


#' Escape special regex characters (Internal)
#' @param x Character string to escape
#' @return Escaped string safe for regex
#' @keywords internal
.escape_regex_fixed <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub(".", "\\.", x, fixed = TRUE)
  x <- gsub("|", "\\|", x, fixed = TRUE)
  x <- gsub("(", "\\(", x, fixed = TRUE)
  x <- gsub(")", "\\)", x, fixed = TRUE)
  x <- gsub("[", "\\[", x, fixed = TRUE)
  x <- gsub("]", "\\]", x, fixed = TRUE)
  x <- gsub("{", "\\{", x, fixed = TRUE)
  x <- gsub("}", "\\}", x, fixed = TRUE)
  x <- gsub("^", "\\^", x, fixed = TRUE)
  x <- gsub("$", "\\$", x, fixed = TRUE)
  x <- gsub("*", "\\*", x, fixed = TRUE)
  x <- gsub("+", "\\+", x, fixed = TRUE)
  x <- gsub("?", "\\?", x, fixed = TRUE)
  return(x)
}
