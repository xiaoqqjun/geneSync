# geneSync <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/xiaoqqjun/geneSync/workflows/R-CMD-check/badge.svg)](https://github.com/xiaoqqjun/geneSync/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/geneSync)](https://CRAN.R-project.org/package=geneSync)
<!-- badges: end -->

## Overview

**geneSync** provides fast and reliable gene symbol-to-ID conversion with comprehensive synonym support and cross-species ortholog mapping. It solves common problems in bioinformatics: inconsistent gene naming across different data sources and cross-species gene conversion.

### Key Features

- 🔄 **Bidirectional conversion**: Symbol → ID and ID → Symbol
- 📚 **Synonym support**: Matches historical gene names automatically
- 🧬 **Authority symbols**: Proper handling of mitochondrial genes (MT-ATP6 vs ATP6)
- 🐭 **Multi-species**: Human, Mouse, and Rat databases included
- 🔀 **Ortholog mapping**: Cross-species gene conversion (NEW in v2.3.0)
- 🧫 **Seurat integration**: Direct gene sync for Seurat objects (NEW in v2.3.0)
- 📊 **Detailed statistics**: Comprehensive conversion summary and logging
- 📦 **Self-contained**: No external API calls, works offline
- ⚡ **Fast**: Vectorized operations for rapid lookup

## Installation

```r
# Install from GitHub
devtools::install_github("xiaoqqjun/geneSync")
```

## Quick Start

### Gene Symbol Conversion

```r
library(geneSync)

# Convert gene symbols to IDs
result <- gene_convert(c("GAPDH", "TP53", "BRCA1", "ATP6"), species = "homo")
result
#>   input symbol symbol_authority gene_id      ensembl_id chromosome      match_type final_symbol
#> 1 GAPDH  GAPDH            GAPDH    2597 ENSG00000111640         12 exact_authority        GAPDH
#> 2  TP53   TP53             TP53    7157 ENSG00000141510         17 exact_authority         TP53
#> 3 BRCA1  BRCA1            BRCA1     672 ENSG00000012048         17 exact_authority        BRCA1
#> 4  ATP6   ATP6          MT-ATP6    4508            <NA>         MT    exact_symbol      MT-ATP6

# View conversion summary
convert_summary(result)

# Save log to file
convert_summary(result, save_log = TRUE, log_file = "conversion_log.txt")
```

### Cross-Species Ortholog Conversion (NEW!)

```r
# Human to Mouse
gene_ortholog(c("GAPDH", "TP53", "BRCA1"), from = "homo", to = "mmu")
#>   input homo_symbol homo_symbol_authority homo_gene_id mmu_symbol mmu_symbol_authority mmu_gene_id match_type
#> 1 GAPDH       GAPDH                 GAPDH         2597      Gapdh                Gapdh       14433      found
#> 2  TP53        TP53                  TP53         7157      Trp53                Trp53       22059      found
#> 3 BRCA1       BRCA1                 BRCA1          672      Brca1                Brca1       12189      found

# Mouse to Human
gene_ortholog(c("Gapdh", "Trp53"), from = "mmu", to = "homo")

# Human to Rat
gene_ortholog(c("GAPDH", "TP53"), from = "homo", to = "rat")

# Mouse to Rat
gene_ortholog(c("Gapdh", "Trp53"), from = "mmu", to = "rat")

# Batch conversion with statistics
result <- batch_ortholog(genes, from = "homo", to = "mmu", verbose = TRUE)
```

### Seurat Object Integration (NEW!)

```r
library(Seurat)
library(geneSync)

# Within-species gene conversion for human scRNA-seq data
obj <- gene_sync_add_to_obj(obj, species = "homo", verbose = TRUE)

# View conversion results
sync_table <- get_gene_sync_table(obj)
head(sync_table)

# Get only kept genes
kept_genes <- get_gene_sync_table(obj, filter = "yes")

# Get removed genes
removed_genes <- get_gene_sync_table(obj, filter = "no")

# Print detailed summary statistics
print_gene_sync_summary(obj)

# Save summary to file
print_gene_sync_summary(obj, save_log = TRUE, log_file = "sync_summary.txt")

# Export sync table to file
export_gene_sync_table(obj, file_path = "gene_sync_table.tsv")

# Export only kept genes
export_gene_sync_table(obj, file_path = "kept_genes.tsv", filter = "yes")
```

#### Keep Status Explanation

| Status | Description |
|--------|-------------|
| `yes` | Gene is reliable and retained in the object |
| `no` | Gene name is unreliable and filtered out |
| `others` | Gene needs further validation |

#### Cross-species conversion for Seurat objects

```r
# Human to Mouse ortholog mapping
obj <- gene_sync_add_to_obj(obj, species = "homo", ortholog_to = "mmu")
```

## Functions

### Gene Conversion

| Function | Description |
|----------|-------------|
| `gene_convert()` | Convert gene symbols ↔ IDs |
| `gene_db_info()` | Display gene database information |
| `convert_summary()` | Generate conversion statistics (with logging) |
| `export_unmatched()` | Export unmatched genes |
| `export_converted()` | Export synonym-converted genes |
| `get_mito_genes()` | Get mitochondrial gene list |
| `get_chr_genes()` | Get genes by chromosome |

### Ortholog Conversion

| Function | Description |
|----------|-------------|
| `gene_ortholog()` | Convert genes between species |
| `batch_ortholog()` | Batch conversion with statistics |
| `ortholog_db_info()` | Display ortholog database information |

### Seurat Object Integration

| Function | Description |
|----------|-------------|
| `gene_sync_add_to_obj()` | Sync gene names in Seurat objects |
| `get_gene_sync_table()` | Extract gene sync table from Seurat object |
| `print_gene_sync_summary()` | Display gene sync statistics |
| `export_gene_sync_table()` | Export gene sync table to file |

## Database Information

### Gene Databases

```r
gene_db_info()
#> ==================================================
#>               geneSync Database Info
#> ==================================================
#> Data Source    : NCBI Gene
#> Update Date    : 2026-02-05
#>
#> [ Homo sapiens (Human) ]
#>   Total genes        : 193,859
#>   Protein-coding     : 20,xxx
#>   With authority sym : 44,938 (23.2%)
#>   Mitochondrial      : 37
```

### Ortholog Databases

```r
ortholog_db_info()
#> ==================================================
#>            geneSync Ortholog Database
#> ==================================================
#> Data Source    : NCBI gene_orthologs
#> Update Date    : 2026-02-05
#>
#> Human - Mouse  : 16,822 pairs
#> Human - Rat    : 16,753 pairs
#> Mouse - Rat    : 16,494 pairs (bridged via human)
#>
#> - All mappings are one-to-one orthologs
```

## Output Columns

### gene_convert()

| Column | Description |
|--------|-------------|
| input | Original input gene name |
| symbol | NCBI gene symbol |
| symbol_authority | Official symbol (HGNC/MGI/RGD), NA if not available |
| gene_id | Entrez Gene ID |
| ensembl_id | Ensembl gene ID |
| chromosome | Chromosome location |
| match_type | How the match was made |
| final_symbol | Recommended symbol: authority > symbol > input |

### get_gene_sync_table() (Seurat)

| Column | Description |
|--------|-------------|
| orig.feature | Original feature name from input |
| new.feature | Updated feature name (final_symbol or ortholog symbol) |
| symbol | Official gene symbol |
| symbol_authority | Authoritative symbol (e.g., HGNC approved) |
| gene_id | NCBI Gene ID |
| ensembl_id | Ensembl ID |
| chromosome | Chromosome location |
| match_type | Type of match (exact_authority, exact_symbol, synonym, not_found, etc.) |
| final_symbol | Final symbol used |
| keep | Whether the gene was kept (yes/no/others) |
| duplication | Whether there was duplication (yes/no) |
| n_final | Number of inputs mapping to this final_symbol |

### gene_ortholog()

| Column | Description |
|--------|-------------|
| input | Original input gene |
| [from]_symbol | Source species symbol |
| [from]_symbol_authority | Source species authority symbol |
| [from]_gene_id | Source species Entrez ID |
| [to]_symbol | Target species symbol |
| [to]_symbol_authority | Target species authority symbol |
| [to]_gene_id | Target species Entrez ID |
| match_type | "found" or "not_found" |

## Match Types

| Type | Description |
|------|-------------|
| `exact_authority` | Matched official authority symbol |
| `exact_symbol` | Matched NCBI symbol |
| `synonym` | Matched via historical/alternative name |
| `id_match` | Matched by Entrez Gene ID |
| `not_found` | No match found |
| `invalid_id` | Invalid numeric ID |

## Citation

If you use geneSync in your research, please cite:

```
Feng Z (2026). geneSync: Gene Symbol and ID Conversion with Synonym Support.
R package version 2.3.0. https://github.com/xiaoqqjun/geneSync
```

## Author

**Zhijun Feng** (Leone)

- Email: xiaoqqjun@sina.com
- ORCID: [0000-0003-1813-1669](https://orcid.org/0000-0003-1813-1669)
- GitHub: [@xiaoqqjun](https://github.com/xiaoqqjun)
- WeChat: 博士后的小酒馆

## License

MIT License
