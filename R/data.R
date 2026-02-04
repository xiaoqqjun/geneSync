#' Human Gene Database
#'
#' A dataset containing gene information for Homo sapiens (human) from NCBI Gene.
#'
#' @format A data frame with 8 columns:
#' \describe{
#'   \item{gene_id}{NCBI Entrez Gene ID}
#'   \item{symbol}{NCBI gene symbol}
#'   \item{symbol_authority}{Official symbol from HGNC (Human Gene Nomenclature Committee),
#'     or NA if not available. For mitochondrial genes, this contains the standard
#'     MT- prefixed names (e.g., MT-ATP6)}
#'   \item{synonyms}{Alternative gene names separated by "|", or "-" if none}
#'   \item{ensembl_id}{Ensembl gene ID, or NA if not available}
#'   \item{chromosome}{Chromosome location (e.g., "1", "X", "MT")}
#'   \item{type}{Gene type (e.g., "protein-coding", "pseudo", "ncRNA")}
#'   \item{description}{Gene description/full name}
#' }
#'
#' @source NCBI Gene \url{https://www.ncbi.nlm.nih.gov/gene}
#' @seealso \code{\link{gene_convert}}, \code{\link{gene_db_info}}
"gene_db_homo"


#' Mouse Gene Database
#'
#' A dataset containing gene information for Mus musculus (mouse) from NCBI Gene.
#'
#' @format A data frame with 8 columns:
#' \describe{
#'   \item{gene_id}{NCBI Entrez Gene ID}
#'   \item{symbol}{NCBI gene symbol}
#'   \item{symbol_authority}{Official symbol from MGI (Mouse Genome Informatics),
#'     or NA if not available. For mitochondrial genes, this contains the standard
#'     mt- prefixed names (e.g., mt-Atp6)}
#'   \item{synonyms}{Alternative gene names separated by "|", or "-" if none}
#'   \item{ensembl_id}{Ensembl gene ID, or NA if not available}
#'   \item{chromosome}{Chromosome location (e.g., "1", "X", "MT")}
#'   \item{type}{Gene type (e.g., "protein-coding", "pseudo", "ncRNA")}
#'   \item{description}{Gene description/full name}
#' }
#'
#' @source NCBI Gene \url{https://www.ncbi.nlm.nih.gov/gene}
#' @seealso \code{\link{gene_convert}}, \code{\link{gene_db_info}}
"gene_db_mmu"


#' Rat Gene Database
#'
#' A dataset containing gene information for Rattus norvegicus (rat) from NCBI Gene.
#'
#' @format A data frame with 8 columns:
#' \describe{
#'   \item{gene_id}{NCBI Entrez Gene ID}
#'   \item{symbol}{NCBI gene symbol}
#'   \item{symbol_authority}{Official symbol from RGD (Rat Genome Database),
#'     or NA if not available. For mitochondrial genes, this contains the standard
#'     Mt- prefixed names}
#'   \item{synonyms}{Alternative gene names separated by "|", or "-" if none}
#'   \item{ensembl_id}{Ensembl gene ID, or NA if not available}
#'   \item{chromosome}{Chromosome location (e.g., "1", "X", "MT")}
#'   \item{type}{Gene type (e.g., "protein-coding", "pseudo", "ncRNA")}
#'   \item{description}{Gene description/full name}
#' }
#'
#' @source NCBI Gene \url{https://www.ncbi.nlm.nih.gov/gene}
#' @seealso \code{\link{gene_convert}}, \code{\link{gene_db_info}}
"gene_db_rat"
