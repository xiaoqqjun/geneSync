#' Conversion Summary Statistics
#'
#' Generate detailed summary statistics for gene conversion results.
#'
#' @param result Data frame returned by \code{gene_convert()}.
#' @param save_log Logical. If TRUE, save the summary to a text file. Default is FALSE.
#' @param log_file Character. Path to save the log file. Default is "geneSync_log.txt".
#'
#' @return Invisibly returns a list containing all statistics.
#'
#' @details
#' The summary includes:
#' \itemize{
#'   \item Total input genes
#'   \item Match type breakdown (exact_authority, exact_symbol, synonym, not_found)
#'   \item Consistent genes: input = symbol = symbol_authority (excluding NA)
#'   \item Converted genes: input differs from symbol but symbol = symbol_authority (excluding NA)
#'   \item Unmatched genes: no symbol or symbol_authority found
#' }
#'
#' @examples
#' \dontrun{
#' result <- gene_convert(genes, species = "homo")
#' convert_summary(result)
#'
#' # Save log to file
#' convert_summary(result, save_log = TRUE, log_file = "my_conversion_log.txt")
#' }
#'
#' @export
convert_summary <- function(result, save_log = FALSE, log_file = "geneSync_log.txt") {
  
  # Validate input
  required_cols <- c("input", "symbol", "symbol_authority", "match_type")
  if (!all(required_cols %in% names(result))) {
    stop("Input must be a result from gene_convert(). Missing columns: ",
         paste(setdiff(required_cols, names(result)), collapse = ", "),
         call. = FALSE)
  }
  
  n_total <- nrow(result)
  
  # ========================================
  # Match type statistics
  # ========================================
  match_type_table <- table(result$match_type, useNA = "ifany")
  
  n_exact_authority <- sum(result$match_type == "exact_authority", na.rm = TRUE)
  n_exact_symbol <- sum(result$match_type == "exact_symbol", na.rm = TRUE)
  n_synonym <- sum(result$match_type == "synonym", na.rm = TRUE)
  n_not_found <- sum(result$match_type == "not_found", na.rm = TRUE)
  n_id_match <- sum(result$match_type == "id_match", na.rm = TRUE)
  n_invalid_id <- sum(result$match_type == "invalid_id", na.rm = TRUE)
  
  n_matched <- n_total - n_not_found - n_invalid_id
  
  # ========================================
  # Consistency analysis (excluding NA)
  # ========================================
  
  # Subset: has both symbol and symbol_authority
  has_both <- !is.na(result$symbol) & !is.na(result$symbol_authority)
  subset_both <- result[has_both, ]
  
  # Consistent: input = symbol = symbol_authority (case-insensitive)
  if (nrow(subset_both) > 0) {
    consistent_idx <- toupper(subset_both$input) == toupper(subset_both$symbol) &
                      toupper(subset_both$symbol) == toupper(subset_both$symbol_authority)
    n_all_consistent <- sum(consistent_idx)
    
    # Converted: input != symbol, but symbol = symbol_authority
    converted_idx <- toupper(subset_both$input) != toupper(subset_both$symbol) &
                     toupper(subset_both$symbol) == toupper(subset_both$symbol_authority)
    n_converted <- sum(converted_idx)
    
    # input matches symbol_authority but not symbol (e.g., MT-ATP6 input, ATP6 symbol, MT-ATP6 authority)
    auth_match_idx <- toupper(subset_both$input) == toupper(subset_both$symbol_authority) &
                      toupper(subset_both$input) != toupper(subset_both$symbol)
    n_auth_match <- sum(auth_match_idx)
  } else {
    n_all_consistent <- 0
    n_converted <- 0
    n_auth_match <- 0
  }
  
  # Subset: has symbol only (no authority)
  has_symbol_only <- !is.na(result$symbol) & is.na(result$symbol_authority)
  n_symbol_only <- sum(has_symbol_only)
  
  # Unmatched: no symbol and no symbol_authority
  unmatched_idx <- is.na(result$symbol) & is.na(result$symbol_authority)
  n_unmatched <- sum(unmatched_idx)
  
  # ========================================
  # Build log content
  # ========================================
  
  log_lines <- c(
    "========================================================",
    "           geneSync Conversion Summary Report           ",
    "========================================================",
    paste0("Generated: ", Sys.time()),
    paste0("geneSync version: ", as.character(utils::packageVersion("geneSync"))),
    "",
    "--------------------------------------------------------",
    "                   OVERALL STATISTICS                   ",
    "--------------------------------------------------------",
    sprintf("Total input genes        : %s", format(n_total, big.mark = ",")),
    sprintf("Successfully matched     : %s (%.2f%%)", 
            format(n_matched, big.mark = ","), n_matched / n_total * 100),
    sprintf("Not found                : %s (%.2f%%)", 
            format(n_not_found, big.mark = ","), n_not_found / n_total * 100),
    "",
    "--------------------------------------------------------",
    "                 MATCH TYPE BREAKDOWN                   ",
    "--------------------------------------------------------",
    sprintf("exact_authority          : %s (%.2f%%)", 
            format(n_exact_authority, big.mark = ","), n_exact_authority / n_total * 100),
    sprintf("exact_symbol             : %s (%.2f%%)", 
            format(n_exact_symbol, big.mark = ","), n_exact_symbol / n_total * 100),
    sprintf("synonym                  : %s (%.2f%%)", 
            format(n_synonym, big.mark = ","), n_synonym / n_total * 100),
    sprintf("id_match                 : %s (%.2f%%)", 
            format(n_id_match, big.mark = ","), n_id_match / n_total * 100),
    sprintf("not_found                : %s (%.2f%%)", 
            format(n_not_found, big.mark = ","), n_not_found / n_total * 100),
    sprintf("invalid_id               : %s (%.2f%%)", 
            format(n_invalid_id, big.mark = ","), n_invalid_id / n_total * 100),
    "",
    "--------------------------------------------------------",
    "               CONSISTENCY ANALYSIS                     ",
    "--------------------------------------------------------",
    sprintf("All consistent           : %s", format(n_all_consistent, big.mark = ",")),
    "  (input = symbol = symbol_authority, excluding NA)",
    "",
    sprintf("Converted via synonym    : %s", format(n_converted, big.mark = ",")),
    "  (input != symbol, but symbol = symbol_authority)",
    "",
    sprintf("Authority matched        : %s", format(n_auth_match, big.mark = ",")),
    "  (input = symbol_authority != symbol, e.g., MT-genes)",
    "",
    sprintf("Symbol only (no auth)    : %s", format(n_symbol_only, big.mark = ",")),
    "  (matched symbol but symbol_authority is NA)",
    "",
    sprintf("Unmatched (input only)   : %s", format(n_unmatched, big.mark = ",")),
    "  (no symbol or symbol_authority found)",
    "",
    "--------------------------------------------------------",
    "                      NOTES                             ",
    "--------------------------------------------------------",
    "- final_symbol column: symbol_authority > symbol > input",
    "- Use match_type column to filter specific categories",
    "- Unmatched genes are preserved with original input name",
    "",
    "========================================================",
    "Data source: NCBI Gene | https://github.com/xiaoqqjun/geneSync",
    "========================================================"
  )
  
  # Print to console
  cat(paste(log_lines, collapse = "\n"), "\n")
  
  # Save to file if requested
  if (save_log) {
    writeLines(log_lines, con = log_file)
    message("\nLog saved to: ", log_file)
  }
  
  # Return statistics invisibly
  stats <- list(
    total = n_total,
    matched = n_matched,
    not_found = n_not_found,
    match_type = list(
      exact_authority = n_exact_authority,
      exact_symbol = n_exact_symbol,
      synonym = n_synonym,
      id_match = n_id_match,
      invalid_id = n_invalid_id
    ),
    consistency = list(
      all_consistent = n_all_consistent,
      converted = n_converted,
      auth_match = n_auth_match,
      symbol_only = n_symbol_only,
      unmatched = n_unmatched
    )
  )
  
  invisible(stats)
}


#' Export Unmatched Genes
#'
#' Extract genes that could not be matched from conversion results.
#'
#' @param result Data frame returned by \code{gene_convert()}.
#' @param save_file Logical. If TRUE, save to a text file. Default is FALSE.
#' @param file_path Character. Path to save the file. Default is "unmatched_genes.txt".
#'
#' @return Character vector of unmatched gene names.
#'
#' @examples
#' \dontrun{
#' result <- gene_convert(genes, species = "homo")
#' unmatched <- export_unmatched(result, save_file = TRUE)
#' }
#'
#' @export
export_unmatched <- function(result, save_file = FALSE, file_path = "unmatched_genes.txt") {
  
  unmatched_genes <- result$input[result$match_type == "not_found"]
  
  message("Unmatched genes: ", length(unmatched_genes))
  
  if (save_file && length(unmatched_genes) > 0) {
    writeLines(unmatched_genes, con = file_path)
    message("Saved to: ", file_path)
  }
  
  return(unmatched_genes)
}


#' Export Converted Genes (Synonym Matches)
#'
#' Extract genes that were matched via synonyms (i.e., input name differs from official symbol).
#'
#' @param result Data frame returned by \code{gene_convert()}.
#' @param save_file Logical. If TRUE, save to a text file. Default is FALSE.
#' @param file_path Character. Path to save the file. Default is "converted_genes.txt".
#'
#' @return Data frame with input and final_symbol for converted genes.
#'
#' @examples
#' \dontrun{
#' result <- gene_convert(genes, species = "homo")
#' converted <- export_converted(result, save_file = TRUE)
#' }
#'
#' @export
export_converted <- function(result, save_file = FALSE, file_path = "converted_genes.txt") {
  
  # Genes where input differs from final_symbol
  converted_idx <- result$match_type != "not_found" &
                   toupper(result$input) != toupper(result$final_symbol)
  
  converted_df <- result[converted_idx, c("input", "symbol", "symbol_authority", "final_symbol", "match_type")]
  
  message("Converted genes: ", nrow(converted_df))
  
  if (save_file && nrow(converted_df) > 0) {
    write.table(converted_df, file = file_path, sep = "\t", 
                row.names = FALSE, quote = FALSE)
    message("Saved to: ", file_path)
  }
  
  return(converted_df)
}
