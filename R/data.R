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


#' Human-Mouse Ortholog Mapping
#'
#' A dataset containing one-to-one ortholog mappings between human and mouse genes.
#'
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{homo_gene_id}{Human Entrez Gene ID}
#'   \item{homo_symbol}{Human NCBI gene symbol}
#'   \item{homo_symbol_authority}{Human HGNC official symbol}
#'   \item{mmu_gene_id}{Mouse Entrez Gene ID}
#'   \item{mmu_symbol}{Mouse NCBI gene symbol}
#'   \item{mmu_symbol_authority}{Mouse MGI official symbol}
#' }
#'
#' @source NCBI gene_orthologs \url{https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz}
#' @seealso \code{\link{gene_ortholog}}, \code{\link{ortholog_db_info}}
"ortholog_homo_mmu"


#' Human-Rat Ortholog Mapping
#'
#' A dataset containing one-to-one ortholog mappings between human and rat genes.
#'
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{homo_gene_id}{Human Entrez Gene ID}
#'   \item{homo_symbol}{Human NCBI gene symbol}
#'   \item{homo_symbol_authority}{Human HGNC official symbol}
#'   \item{rat_gene_id}{Rat Entrez Gene ID}
#'   \item{rat_symbol}{Rat NCBI gene symbol}
#'   \item{rat_symbol_authority}{Rat RGD official symbol}
#' }
#'
#' @source NCBI gene_orthologs \url{https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz}
#' @seealso \code{\link{gene_ortholog}}, \code{\link{ortholog_db_info}}
"ortholog_homo_rat"


#' Mouse-Rat Ortholog Mapping
#'
#' A dataset containing one-to-one ortholog mappings between mouse and rat genes,
#' derived through human bridging.
#'
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{mmu_gene_id}{Mouse Entrez Gene ID}
#'   \item{mmu_symbol}{Mouse NCBI gene symbol}
#'   \item{mmu_symbol_authority}{Mouse MGI official symbol}
#'   \item{rat_gene_id}{Rat Entrez Gene ID}
#'   \item{rat_symbol}{Rat NCBI gene symbol}
#'   \item{rat_symbol_authority}{Rat RGD official symbol}
#' }
#'
#' @source NCBI gene_orthologs (bridged through human)
#' @seealso \code{\link{gene_ortholog}}, \code{\link{ortholog_db_info}}
"ortholog_mmu_rat"
