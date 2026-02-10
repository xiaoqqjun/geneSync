# Test loading ortholog file
tryCatch({
  load('ortholog_homo_mmu.rda')
  cat('Ortholog file loaded successfully!\n')
  cat('Object names:', ls(), '\n')
}, error = function(e) {
  cat('Error loading ortholog file:', e$message, '\n')
})
