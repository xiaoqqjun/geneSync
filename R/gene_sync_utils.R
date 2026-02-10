#' Get Gene Sync Table from Seurat Object
#'
#' Extract the gene synchronization table from a Seurat object processed by
#' \code{gene_sync_add_to_obj()}.
#'
#' @param object A Seurat object that has been processed by \code{gene_sync_add_to_obj()}.
#' @param filter Character. Filter the table by keep status. Options are:
#'   \itemize{
#'     \item "all" (default): Return all rows
#'     \item "yes": Return only keep = "yes" rows
#'     \item "no": Return only keep = "no" rows
#'     \item "others": Return only keep = "others" rows
#'   }
#' @param include_stats Logical. If TRUE, includes summary statistics as an attribute.
#'   Default is FALSE.
#'
#' @return A data frame containing the gene synchronization table.
#'
#' @details
#' The returned table contains the following columns:
#' \itemize{
#'   \item orig.feature: Original feature name from input
#'   \item new.feature: Updated feature name (final_symbol or ortholog symbol)
#'   \item symbol: Official gene symbol
#'   \item symbol_authority: Authoritative symbol (e.g., HGNC approved)
#'   \item gene_id: NCBI Gene ID
#'   \item ensembl_id: Ensembl ID
#'   \item chromosome: Chromosome location
#'   \item match_type: Type of match (exact_authority, exact_symbol, synonym, not_found, etc.)
#'   \item final_symbol: Final symbol used
#'   \item keep: Whether the gene was kept (yes/no/others)
#'   \item duplication: Whether there was duplication (yes/no)
#' }
#'
#' For ortholog mode, columns include species-specific information.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(geneSync)
#'
#' # Process a Seurat object
#' obj <- gene_sync_add_to_obj(obj, species = "homo")
#'
#' # Get the full sync table
#' sync_table <- get_gene_sync_table(obj)
#'
#' # Get only kept genes
#' kept_table <- get_gene_sync_table(obj, filter = "yes")
#'
#' # Get genes that were not kept
#' removed_table <- get_gene_sync_table(obj, filter = "no")
#' }
#'
#' @export
get_gene_sync_table <- function(object, filter = c("all", "yes", "no", "others"),
                                 include_stats = FALSE) {

  # Check if gene_sync info exists
  if (is.null(object@misc$gene_sync)) {
    stop("No gene_sync information found in object@misc$gene_sync. ",
         "Please run gene_sync_add_to_obj() first.",
         call. = FALSE)
  }

  filter <- match.arg(filter)

  # Extract the sync table
  sync_table <- object@misc$gene_sync$table

  # Apply filter if requested
  if (filter != "all") {
    if ("keep" %in% names(sync_table)) {
      sync_table <- sync_table[sync_table$keep == filter, , drop = FALSE]
    } else {
      warning("keep column not found in sync table. Returning all rows.")
    }
  }

  # Add summary stats as attribute if requested
  if (include_stats) {
    attr(sync_table, "summary") <- object@misc$gene_sync$summary
    attr(sync_table, "settings") <- object@misc$gene_sync$settings
  }

  return(sync_table)
}


#' Print Gene Sync Summary
#'
#' Display a detailed summary of gene synchronization results from a Seurat
#' object processed by \code{gene_sync_add_to_obj()}.
#'
#' @param object A Seurat object that has been processed by \code{gene_sync_add_to_obj()}.
#' @param save_log Logical. If TRUE, save the summary to a text file. Default is FALSE.
#' @param log_file Character. Path to save the log file. Default is "gene_sync_summary.txt".
#'
#' @return Invisibly returns a list containing all statistics.
#'
#' @details
#' The summary includes:
#' \itemize{
#'   \item Total input features
#'   \item Final features after filtering
#'   \item Match statistics (matched, not_found)
#'   \item Keep status breakdown (yes, no, others)
#'   \item Duplication statistics
#'   \item Match type distribution
#'   \item Processing settings (species, symbol_type, etc.)
#'   \item Timestamp
#' }
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(geneSync)
#'
#' # Process a Seurat object
#' obj <- gene_sync_add_to_obj(obj, species = "homo")
#'
#' # Print summary to console
#' print_gene_sync_summary(obj)
#'
#' # Save summary to file
#' print_gene_sync_summary(obj, save_log = TRUE, log_file = "my_sync_summary.txt")
#' }
#'
#' @export
print_gene_sync_summary <- function(object, save_log = FALSE,
                                     log_file = "gene_sync_summary.txt") {

  # Check if gene_sync info exists
  if (is.null(object@misc$gene_sync)) {
    stop("No gene_sync information found in object@misc$gene_sync. ",
         "Please run gene_sync_add_to_obj() first.",
         call. = FALSE)
  }

  sync_info <- object@misc$gene_sync
  summary <- sync_info$summary
  settings <- sync_info$settings

  # Build log content
  log_lines <- c(
    "========================================================",
    "         geneSync for Seurat Object Summary            ",
    "========================================================",
    paste0("Generated: ", Sys.time()),
    paste0("geneSync version: ", as.character(utils::packageVersion("geneSync"))),
    "",
    "--------------------------------------------------------",
    "                    PROCESSING INFO                     ",
    "--------------------------------------------------------"
  )

  # Add sync type and species info
  if (sync_info$type == "convert") {
    log_lines <- c(log_lines,
      sprintf("Sync Type                : Gene Conversion"),
      sprintf("Species                  : %s", settings$species),
      sprintf("Symbol Type              : %s", settings$symbol_type),
      sprintf("Update Features          : %s", settings$update_features)
    )
  } else if (sync_info$type == "ortholog") {
    log_lines <- c(log_lines,
      sprintf("Sync Type                : Ortholog Mapping"),
      sprintf("From Species             : %s", settings$from_species),
      sprintf("To Species               : %s", settings$to_species),
      sprintf("Symbol Type              : %s", settings$symbol_type),
      sprintf("Update Features          : %s", settings$update_features)
    )
  }

  log_lines <- c(log_lines,
    sprintf("Timestamp                : %s", as.character(settings$timestamp)),
    "",
    "--------------------------------------------------------",
    "                   OVERALL STATISTICS                   ",
    "--------------------------------------------------------",
    sprintf("Total Input Features     : %s", format(summary$total_input, big.mark = ",")),
    sprintf("Final Features           : %s", format(summary$final_features, big.mark = ",")),
    sprintf("Features Removed         : %s (%.2f%%)",
            format(summary$total_input - summary$final_features, big.mark = ","),
            (summary$total_input - summary$final_features) / summary$total_input * 100),
    "",
    sprintf("Matched                  : %s (%.2f%%)",
            format(summary$matched, big.mark = ","),
            summary$matched / summary$total_input * 100),
    sprintf("Not Found                : %s (%.2f%%)",
            format(summary$not_found, big.mark = ","),
            summary$not_found / summary$total_input * 100),
    "",
    "--------------------------------------------------------",
    "                     KEEP STATUS                         ",
    "--------------------------------------------------------",
    sprintf("keep = yes               : %s (%.2f%%)",
            format(summary$keep_yes, big.mark = ","),
            summary$keep_yes / summary$total_input * 100),
    sprintf("keep = no                : %s (%.2f%%)",
            format(summary$keep_no, big.mark = ","),
            summary$keep_no / summary$total_input * 100)
  )

  if (summary$keep_others > 0) {
    log_lines <- c(log_lines,
      sprintf("keep = others            : %s (%.2f%%)",
              format(summary$keep_others, big.mark = ","),
              summary$keep_others / summary$total_input * 100)
    )
  }

  log_lines <- c(log_lines,
    "",
    "--------------------------------------------------------",
    "                   DUPLICATION INFO                     ",
    "--------------------------------------------------------",
    sprintf("duplication = yes        : %s",
            format(summary$duplication_yes, big.mark = ","))
  )

  # Add match type breakdown if available
  if (!is.null(summary$match_type_counts)) {
    log_lines <- c(log_lines,
      "",
      "--------------------------------------------------------",
    "                 MATCH TYPE BREAKDOWN                   ",
      "--------------------------------------------------------"
    )

    match_types <- names(summary$match_type_counts)
    for (mt in match_types) {
      count <- summary$match_type_counts[[mt]]
      pct <- count / summary$total_input * 100
      log_lines <- c(log_lines,
        sprintf("%-25s : %s (%.2f%%)",
                paste0("  ", mt),
                format(count, big.mark = ","),
                pct)
      )
    }
  }

  log_lines <- c(log_lines,
    "",
    "--------------------------------------------------------",
    "                      NOTES                             ",
    "--------------------------------------------------------",
    "- keep = yes: Gene is reliable and retained",
    "- keep = no: Gene name is unreliable and filtered out",
    "- keep = others: Gene needs further validation",
    "- duplication = yes: Multiple inputs map to same final_symbol",
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
  invisible(list(
    summary = summary,
    settings = settings,
    sync_type = sync_info$type
  ))
}


#' Export Gene Sync Table
#'
#' Export the gene synchronization table to a file.
#'
#' @param object A Seurat object that has been processed by \code{gene_sync_add_to_obj()}.
#' @param file_path Character. Path to save the file. Default is "gene_sync_table.tsv".
#' @param filter Character. Filter the table by keep status before exporting. Default is "all".
#'
#' @return Invisibly returns the sync table.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(geneSync)
#'
#' # Process a Seurat object
#' obj <- gene_sync_add_to_obj(obj, species = "homo")
#'
#' # Export full table
#' export_gene_sync_table(obj, file_path = "my_sync_table.tsv")
#'
#' # Export only kept genes
#' export_gene_sync_table(obj, file_path = "kept_genes.tsv", filter = "yes")
#' }
#'
#' @export
export_gene_sync_table <- function(object, file_path = "gene_sync_table.tsv",
                                    filter = c("all", "yes", "no", "others")) {

  sync_table <- get_gene_sync_table(object, filter = filter)

  write.table(sync_table, file = file_path, sep = "\t",
              row.names = FALSE, quote = FALSE)

  message("Gene sync table exported to: ", file_path)
  message("  Rows: ", nrow(sync_table))

  invisible(sync_table)
}
