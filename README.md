# geneSync <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/xiaoqqjun/geneSync/workflows/R-CMD-check/badge.svg)](https://github.com/xiaoqqjun/geneSync/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/geneSync)](https://CRAN.R-project.org/package=geneSync)
<!-- badges: end -->

## Overview

**geneSync** provides fast and reliable gene symbol-to-ID conversion with comprehensive synonym support. It solves a common problem in bioinformatics: inconsistent gene naming across different data sources and platforms.

### Key Features

- 🔄 **Bidirectional conversion**: Symbol → ID and ID → Symbol
- 📚 **Synonym support**: Matches historical gene names automatically
- 🧬 **Authority symbols**: Proper handling of mitochondrial genes (MT-ATP6 vs ATP6)
- 🐭 **Multi-species**: Human, Mouse, and Rat databases included
- 📊 **Detailed statistics**: Comprehensive conversion summary and logging
- 📦 **Self-contained**: No external API calls, works offline
- ⚡ **Fast**: Vectorized operations for rapid lookup

## Installation

```r
# Install from GitHub
devtools::install_github("xiaoqqjun/geneSync")
```

## Quick Start

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

## Output Columns

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

## Match Types

| Type | Description |
|------|-------------|
| `exact_authority` | Matched official authority symbol |
| `exact_symbol` | Matched NCBI symbol |
| `synonym` | Matched via historical/alternative name |
| `id_match` | Matched by Entrez Gene ID |
| `not_found` | No match found |
| `invalid_id` | Invalid numeric ID |

## Conversion Summary

```r
convert_summary(result)
#> ========================================================
#>            geneSync Conversion Summary Report           
#> ========================================================
#> 
#> OVERALL STATISTICS
#> --------------------------------------------------------
#> Total input genes        : 20,000
#> Successfully matched     : 18,500 (92.50%)
#> Not found                : 1,500 (7.50%)
#> 
#> MATCH TYPE BREAKDOWN
#> --------------------------------------------------------
#> exact_authority          : 15,000 (75.00%)
#> exact_symbol             : 3,000 (15.00%)
#> synonym                  : 500 (2.50%)
#> ...
#> 
#> CONSISTENCY ANALYSIS
#> --------------------------------------------------------
#> All consistent           : 14,500
#>   (input = symbol = symbol_authority)
#> Converted via synonym    : 500
#>   (input != symbol, but symbol = symbol_authority)
#> ...
```

## Export Functions

```r
# Export unmatched genes for manual review
unmatched <- export_unmatched(result, save_file = TRUE, file_path = "unmatched.txt")

# Export genes that were converted (synonym matches)
converted <- export_converted(result, save_file = TRUE, file_path = "converted.txt")
```

## Functions

| Function | Description |
|----------|-------------|
| `gene_convert()` | Convert gene symbols ↔ IDs |
| `gene_db_info()` | Display database information |
| `convert_summary()` | Generate conversion statistics (with logging) |
| `export_unmatched()` | Export unmatched genes |
| `export_converted()` | Export synonym-converted genes |
| `get_mito_genes()` | Get mitochondrial gene list |
| `get_chr_genes()` | Get genes by chromosome |

## Database Information

```r
gene_db_info()
#> ==================================================
#>               geneSync Database Info              
#> ==================================================
#> Data Source    : NCBI Gene
#> Update Date    : 2026-02-04
#> 
#> [ Homo sapiens (Human) ]
#>   Authority          : HGNC
#>   Total genes        : 193,859
#>   With authority sym : 44,938 (23.2%)
#>   Mitochondrial      : 37
```

## Citation

If you use geneSync in your research, please cite:

```
Feng Z (2026). geneSync: Gene Symbol and ID Conversion with Synonym Support.
R package version 2.2.0. https://github.com/xiaoqqjun/geneSync
```

## Author

**Zhijun Feng** (Leone)

- Email: xiaoqqjun@sina.com
- ORCID: [0000-0003-1813-1669](https://orcid.org/0000-0003-1813-1669)
- GitHub: [@xiaoqqjun](https://github.com/xiaoqqjun)
- WeChat: 博士后的小酒馆

## License

MIT License
