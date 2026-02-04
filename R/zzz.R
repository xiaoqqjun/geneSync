.onAttach <- function(libname, pkgname) {
  
  version <- as.character(utils::packageVersion("geneSync"))
  
  msg <- paste0(
    "\n",
    "========================================\n",
    "        geneSync v", version, "\n",
    "========================================\n",
    "Data source : NCBI Gene\n",
    "Update date : 2026-02-04\n",
    "Species     : Human / Mouse / Rat\n",
    "----------------------------------------\n",
    "GitHub: https://github.com/xiaoqqjun/geneSync\n",
    "========================================\n",
    "\n",
    "Citation:\n",
    "  Feng Z (2026). geneSync: Gene Symbol and ID\n",
    "  Conversion with Synonym Support. R package\n",
    "  version ", version, ".\n",
    "  https://github.com/xiaoqqjun/geneSync\n",
    "\n",
    "Quick start:\n",
    "  gene_convert(genes, species = \"homo\")\n",
    "  gene_db_info()\n",
    "  convert_summary(result, save_log = TRUE)\n"
  )
  
  packageStartupMessage(msg)
}
