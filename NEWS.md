# geneSync News

## Version 3.2.1 (2026-04-14)

### Bug Fixes
- Fixed Seurat V5 compatibility issues:
  - Corrected `GetAssayData()` to use `layer` parameter instead of deprecated `slot` parameter for Seurat V5 objects
  - Improved Seurat version detection logic to properly distinguish between V4 and V5
  - Fixed "subscript out of bounds" error when processing Seurat V5 objects
  - Fixed "Layer 'data' is empty" warning by implementing proper layer detection
  - Enhanced error handling for layer access in Seurat V5

### Technical Changes
- Updated version detection to check for `@layers` structure in addition to layer existence
- Improved subset operation for Seurat V5 layers to ensure gene names exist in layer rownames
- Better fallback logic when primary layers are unavailable

## Version 3.2.0 (2026-02-25)

### New Features
- Added `from` and `to` parameters for flexible input/output format support
- Added `symbol_type` parameter to choose between different symbol columns
- Added `keep` column (yes/no/others) for gene quality classification
- Added `duplication` column (yes/no) for duplicate gene detection
- Implemented detection rate-based duplicate gene resolution
- Added `update_features` parameter to preview changes without modifying object

### Improvements
- Better handling of duplicate gene mappings
- Enhanced statistics and logging
- Support for ortholog mapping between species
