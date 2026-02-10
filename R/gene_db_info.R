#' Display Gene Database Information
#'
#' Show version information and statistics for the built-in gene databases.
#'
#' @param species Optional. If specified, show info for only that species.
#'   Must be one of "homo", "mmu", or "rat". If NULL (default), shows all.
#'
#' @return Invisibly returns a list containing database statistics.
#'
#' @examples
#' # Show all database info
#' gene_db_info()
#'
#' # Show only human database info
#' gene_db_info("homo")
#'
#' @export
gene_db_info <- function(species = NULL) {

  # Database metadata
  db_info <- list(
    homo = list(
      name = "Homo sapiens (Human)",
      authority = "HGNC",
      data = gene_db_homo,
      code = "homo"
    ),
    mmu = list(
      name = "Mus musculus (Mouse)",
      authority = "MGI",
      data = gene_db_mmu,
      code = "mmu"
    ),
    rat = list(
      name = "Rattus norvegicus (Rat)",
      authority = "RGD",
      data = gene_db_rat,
      code = "rat"
    )
  )

  # Filter by species if specified
  if (!is.null(species)) {
    species <- match.arg(species, c("homo", "mmu", "rat"))
    db_info <- db_info[species]
  }

  # Print header
  cat("\n")
  cat("==================================================\n")
  cat("              geneSync Database Info              \n")
  cat("==================================================\n\n")
  cat("Data Source    : NCBI Gene\n")
  cat("Update Date    : 2026-02-04\n")
  cat("Package Version:", as.character(utils::packageVersion("geneSync")), "\n\n")

  # Collect stats
  stats <- list()

  for (sp in names(db_info)) {
    info <- db_info[[sp]]
    db <- info$data

    n_genes <- nrow(db)
    n_with_authority <- sum(!is.na(db$symbol_authority))
    n_with_ensembl <- sum(!is.na(db$ensembl_id))
    n_with_synonyms <- sum(db$synonyms != "-" & !is.na(db$synonyms))
    n_with_chr <- sum(!is.na(db$chromosome))

    # Count by gene type
    type_counts <- table(db$type)
    n_protein_coding <- ifelse("protein-coding" %in% names(type_counts),
                                type_counts["protein-coding"], 0)

    # Count mitochondrial genes
    n_mito <- sum(db$chromosome == "MT", na.rm = TRUE)

    stats[[sp]] <- list(
      total = n_genes,
      protein_coding = n_protein_coding,
      with_authority = n_with_authority,
      with_ensembl = n_with_ensembl,
      with_synonyms = n_with_synonyms,
      mito_genes = n_mito
    )

    cat(sprintf("[ %s ]\n", info$name))
    cat(sprintf("  Species code       : %s\n", info$code))
    cat(sprintf("  Authority          : %s\n", info$authority))
    cat(sprintf("  Total genes        : %s\n", format(n_genes, big.mark = ",")))
    cat(sprintf("  Protein-coding     : %s\n", format(n_protein_coding, big.mark = ",")))
    cat(sprintf("  With authority sym : %s (%.1f%%)\n",
                format(n_with_authority, big.mark = ","),
                n_with_authority / n_genes * 100))
    cat(sprintf("  With Ensembl ID    : %s (%.1f%%)\n",
                format(n_with_ensembl, big.mark = ","),
                n_with_ensembl / n_genes * 100))
    cat(sprintf("  With Synonyms      : %s (%.1f%%)\n",
                format(n_with_synonyms, big.mark = ","),
                n_with_synonyms / n_genes * 100))
    cat(sprintf("  Mitochondrial      : %s\n", format(n_mito, big.mark = ",")))
    cat("\n")
  }

  cat("--------------------------------------------------\n")
  cat("Usage: gene_convert(genes, species = \"homo\")\n")
  cat("==================================================\n\n")

  invisible(stats)
}


#' Get Mitochondrial Genes
#'
#' Retrieve a list of mitochondrial genes for a specified species.
#'
#' @param species Species identifier. Must be one of "homo", "mmu", or "rat".
#' @param return_authority Logical. If TRUE, returns symbol_authority column.
#'   If FALSE (default), returns symbol column.
#'
#' @return Data frame with symbol and symbol_authority for mitochondrial genes.
#'
#' @examples
#' # Get human mitochondrial genes
#' get_mito_genes("homo")
#'
#' @export
get_mito_genes <- function(species = c("homo", "mmu", "rat")) {

  species <- match.arg(species)

  gene_db <- switch(species,
    "homo" = gene_db_homo,
    "mmu"  = gene_db_mmu,
    "rat"  = gene_db_rat
  )

  mito <- gene_db[gene_db$chromosome == "MT" & !is.na(gene_db$chromosome), 
                  c("gene_id", "symbol", "symbol_authority")]

  return(mito)
}


#' Get Genes by Chromosome
#'
#' Retrieve a list of genes on a specified chromosome.
#'
#' @param chr Chromosome identifier (e.g., "1", "X", "MT").
#' @param species Species identifier. Must be one of "homo", "mmu", or "rat".
#'
#' @return Data frame with gene_id, symbol, and symbol_authority.
#'
#' @examples
#' # Get X chromosome genes
#' x_genes <- get_chr_genes("X", species = "homo")
#'
#' # Get mitochondrial genes
#' mt_genes <- get_chr_genes("MT", species = "homo")
#'
#' @export
get_chr_genes <- function(chr, species = c("homo", "mmu", "rat")) {

  species <- match.arg(species)

  gene_db <- switch(species,
    "homo" = gene_db_homo,
    "mmu"  = gene_db_mmu,
    "rat"  = gene_db_rat
  )

  chr_genes <- gene_db[gene_db$chromosome == chr & !is.na(gene_db$chromosome), 
                       c("gene_id", "symbol", "symbol_authority")]

  if (nrow(chr_genes) == 0) {
    warning("No genes found on chromosome: ", chr)
  }

  return(chr_genes)
}
