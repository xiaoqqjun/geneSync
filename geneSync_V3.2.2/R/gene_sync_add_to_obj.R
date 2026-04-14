# gene_sync_add_to_obj 函数 v3.2.2
# 主要更改：
# 1. 添加 from/to 参数，支持多种输入格式
# 2. 添加 update_features 参数（默认 TRUE）
# 3. 添加 keep 列（yes/no/others）
# 4. 添加 duplication 列（yes/no）
# 5. 根据检测率处理重复基因
# 6. 兼容 Seurat V4 和 V5（V5 不支持 drop 参数）
# 7. 修复 Seurat V5 兼容性问题（使用 layer 替代 slot）
# 8. 改进版本检测逻辑，正确处理 Seurat V5 的 layers
# 9. **v3.2.2 修复**:
#    - 修复 Seurat V5: features 必须使用 LogMap() 构造，不能是 factor
#    - 修复 Seurat V5: 处理 layer subset 后变为空矩阵的问题
#    - 修复 Seurat V4: 确保 meta.features 与 data slots 一致

#' Sync Gene Names in Seurat Object
#'
#' Synchronize gene names in a Seurat object to standard symbols.
#' Supports multiple input formats including gene symbols, Ensembl IDs, and Entrez IDs.
#'
#' @param object A Seurat object
#' @param species Source species. Must be one of "homo" (human), "mmu" (mouse), or "rat".
#'   Default is "homo".
#' @param from Input format of gene names in the Seurat object. Must be one of "symbol",
#'   "ensembl_id", or "gene_id". Default is "symbol".
#' @param to Output format for gene names. Must be one of "symbol", "ensembl_id", or "gene_id".
#'   Default is "symbol" (standardizes gene symbols).
#' @param symbol_type Which symbol column to use for final gene names. Must be one of
#'   "final_symbol" (recommended, uses authority symbol if available),
#'   "symbol_authority" (only authority symbols), or "symbol" (NCBI symbols).
#'   Default is "final_symbol".
#' @param ortholog_to Target species for ortholog mapping. If NULL, performs within-species
#'   conversion. Must be one of "homo", "mmu", or "rat" if specified.
#' @param update_features Whether to update the Seurat object features. If TRUE,
#'   removes unmatched genes and renames features. If FALSE, only stores sync info.
#'   Default is TRUE.
#' @param verbose Print progress messages. Default is TRUE.
#'
#' @return A Seurat object with updated gene names and sync information stored in
#'   \code{object@misc$gene_sync}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Convert gene names using gene_convert or gene_ortholog
#'   \item Classify genes based on match quality (keep column)
#'   \item Handle duplicate gene mappings
#'   \item Update Seurat object features if requested
#' }
#'
#' The \code{keep} column in the sync table indicates:
#' \itemize{
#'   \item "yes": Reliable gene, should be kept
#'   \item "no": Unreliable gene name, should be filtered
#'   \item "others": Needs further validation
#' }
#'
#' @examples
#' # Basic usage with gene symbols
#' obj <- gene_sync_add_to_obj(obj, species = "homo")
#'
#' # Convert from Ensembl IDs to symbols
#' obj <- gene_sync_add_to_obj(obj, species = "homo", from = "ensembl_id", to = "symbol")
#'
#' # Convert from Entrez IDs to symbols
#' obj <- gene_sync_add_to_obj(obj, species = "homo", from = "gene_id", to = "symbol")
#'
#' # Cross-species ortholog mapping (mouse to human)
#' obj <- gene_sync_add_to_obj(obj, species = "mmu", ortholog_to = "homo")
#'
#' # Preview without updating features
#' obj <- gene_sync_add_to_obj(obj, species = "homo", update_features = FALSE)
#'
#' @export
gene_sync_add_to_obj <- function(object,
                                  species = c("homo", "mmu", "rat"),
                                  from = c("symbol", "ensembl_id", "gene_id"),
                                  to = c("symbol", "ensembl_id", "gene_id"),
                                  symbol_type = c("final_symbol", "symbol_authority", "symbol"),
                                  ortholog_to = NULL,
                                  update_features = TRUE,
                                  verbose = TRUE) {

  species <- match.arg(species)
  from <- match.arg(from)
  to <- match.arg(to)
  symbol_type <- match.arg(symbol_type)

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (verbose) {
    cat("=== Gene Sync for Seurat Object ===\n")
    cat("Species:", species)
    if (!is.null(ortholog_to)) {
      cat(" ->", ortholog_to)
    }
    cat("\n")
    cat("Input format:", from, "\n")
    cat("Output format:", to, "\n")
    cat("Symbol type:", symbol_type, "\n")
    cat("Update features:", update_features, "\n")
    cat("Total features:", nrow(object), "\n\n")
  }

  # Get original feature names
  orig_features <- rownames(object)
  n_features <- length(orig_features)

  # Perform gene conversion or ortholog mapping
  if (verbose) cat("Processing gene names...\n")
  start_time <- Sys.time()

  if (is.null(ortholog_to)) {
    # Within-species conversion using gene_convert
    conv_result <- gene_convert(orig_features, species = species, from = from, to = to)
    new_features <- conv_result$final_symbol
    sync_type <- "convert"
  } else {
    # Cross-species conversion using gene_ortholog
    ortholog_to <- match.arg(ortholog_to, c("homo", "mmu", "rat"))
    ortho_result <- gene_ortholog(orig_features, from = species, to = ortholog_to)

    # Determine which column to use for new features
    to_symbol_col <- paste0(ortholog_to, "_symbol")
    to_auth_col <- paste0(ortholog_to, "_symbol_authority")

    if (symbol_type == "final_symbol") {
      new_features <- ifelse(!is.na(ortho_result[[to_auth_col]]),
                              ortho_result[[to_auth_col]],
                              ortho_result[[to_symbol_col]])
    } else if (symbol_type == "symbol_authority") {
      new_features <- ortho_result[[to_auth_col]]
    } else {
      new_features <- ortho_result[[to_symbol_col]]
    }

    sync_type <- "ortholog"
  }

  elapsed <- Sys.time() - start_time
  if (verbose) cat("Processing time:", format(elapsed), "\n")

  # Handle unmatched genes (keep original name)
  unmatched <- is.na(new_features)
  new_features[unmatched] <- orig_features[unmatched]

  # ==================== NEW LOGIC FOR keep AND duplication ====================

  if (sync_type == "convert") {
    # Work with a copy for keep/duplication analysis
    genes_table <- conv_result

    # 1. 统计每个 final_symbol 出现的次数
    genes_table <- genes_table %>%
      dplyr::group_by(final_symbol) %>%
      dplyr::mutate(n_final = dplyr::n()) %>%
      dplyr::ungroup()

    # 2. 初始化 keep 为 "others"
    genes_table$keep <- "others"

    # 3. 标记 keep = "yes"
    # Adjust logic based on input type (from parameter)
    if (from == "symbol") {
      genes_table$keep[
        genes_table$n_final == 1 &
          (
            genes_table$match_type == "exact_authority" |
              (
                genes_table$match_type == "synonym" &
                  !is.na(genes_table$symbol_authority) &
                  genes_table$symbol_authority != ""
              )
          )
      ] <- "yes"

      # 4. 标记 keep = "no"
      genes_table$keep[
        genes_table$match_type %in% c("exact_symbol", "not_found") |
          (
            genes_table$match_type == "synonym" &
              (
                is.na(genes_table$symbol_authority) |
                  genes_table$symbol_authority == ""
              )
          )
      ] <- "no"
    } else if (from == "gene_id") {
      # For gene_id input: keep matched IDs
      genes_table$keep[genes_table$match_type == "id_match" & genes_table$n_final == 1] <- "yes"
      genes_table$keep[genes_table$match_type %in% c("not_found", "invalid_id")] <- "no"
    } else if (from == "ensembl_id") {
      # For ensembl_id input: keep matched Ensembl IDs
      genes_table$keep[genes_table$match_type == "ensembl_match" & genes_table$n_final == 1] <- "yes"
      genes_table$keep[genes_table$match_type %in% c("not_found", "invalid_ensembl")] <- "no"
    }

    # 5. 新增 duplication 列
    final_symbol_n <- table(genes_table$final_symbol)
    genes_table$duplication <- ifelse(
      final_symbol_n[genes_table$final_symbol] > 1,
      "yes",
      "no"
    )

    # 6. 处理重复基因 (duplication = "yes" 且 keep = "others")
    dup_rows <- genes_table[genes_table$duplication == "yes" & genes_table$keep == "others", , drop = FALSE]

    if (nrow(dup_rows) > 0) {
      if (verbose) cat("\nDuplicate gene handling:\n")
      if (verbose) cat("  Found", length(unique(dup_rows$final_symbol)), "duplicate final_symbols\n")

      # 获取重复基因的表达数据
      dup_inputs <- unique(dup_rows$input)
      dup_inputs <- intersect(dup_inputs, rownames(object))

      if (length(dup_inputs) > 0) {
        # 使用 GetAssayData 获取表达数据 (V4/V5 兼容)
        # 检测实际对象的类型，而不是包版本
        default_assay <- Seurat::DefaultAssay(object)
        assay_obj <- object[[default_assay]]

        # 改进的 Seurat V5 检测逻辑
        seurat_v5 <- tryCatch({
          # V5 有 @layers 且非空，V4 有 @slots
          has_layers <- !is.null(assay_obj@layers) && length(assay_obj@layers) > 0
          # 额外检查：V5 的 @layers 是一个列表，V4 没有
          if (has_layers && is.list(assay_obj@layers)) {
            TRUE
          } else {
            FALSE
          }
        }, error = function(e) {
          FALSE
        })

        # 获取表达数据
        expr_mat <- NULL
        if (seurat_v5) {
          # Seurat V5: 使用 layer 参数，找到可用的 layer
          layer_names <- names(assay_obj@layers)
          # 在 Seurat V5 中，数据通常在 "data" layer，但也可能在其他 layer
          # 检查哪个 layer 包含数据
          for (target_layer in c("data", "counts", "scale.data")) {
            if (target_layer %in% layer_names) {
              tryCatch({
                layer_data <- assay_obj@layers[[target_layer]]
                if (!is.null(layer_data) && nrow(layer_data) > 0) {
                  # 直接从 layer 获取数据，避免 drop 参数问题
                  expr_mat <- layer_data[dup_inputs, , drop = FALSE]
                  break
                }
              }, error = function(e) {
                # 尝试下一个 layer
              })
            }
          }
          # 如果所有 layer 都失败，尝试使用 GetAssayData
          if (is.null(expr_mat)) {
            tryCatch({
              expr_mat <- Seurat::GetAssayData(object, layer = layer_names[1])
              if (!is.null(dup_inputs) && length(dup_inputs) > 0) {
                expr_mat <- expr_mat[dup_inputs, , drop = FALSE]
              }
            }, error = function(e) {
              if (verbose) cat("  Warning: Could not extract expression data:", conditionMessage(e), "\n")
            })
          }
        } else {
          # Seurat V4: 使用 slot 参数
          tryCatch({
            expr_mat <- Seurat::GetAssayData(object, slot = "data")[dup_inputs, , drop = FALSE]
          }, error = function(e) {
            if (verbose) cat("  Warning: 'data' slot not available, trying 'counts'...\n")
            tryCatch({
              expr_mat <- Seurat::GetAssayData(object, slot = "counts")[dup_inputs, , drop = FALSE]
            }, error = function(e2) {
              if (verbose) cat("  Error: Could not extract expression data:", conditionMessage(e2), "\n")
            })
          })
        }

        mean_exp <- Matrix::rowMeans(expr_mat)
        detect_rate <- Matrix::rowMeans(expr_mat > 0)

        expr_stats <- data.frame(
          input = names(mean_exp),
          mean_expression = mean_exp,
          detection_rate = detect_rate,
          stringsAsFactors = FALSE
        )

        dup_rows2 <- merge(
          dup_rows,
          expr_stats,
          by = "input",
          all.x = TRUE
        )

        # 初始全部标记为 "no"
        dup_rows2$keep <- "no"

        # 对每个 final_symbol，保留检测率最高的
        max_by_group <- ave(
          dup_rows2$detection_rate,
          dup_rows2$final_symbol,
          FUN = max
        )

        dup_rows2$keep[
          dup_rows2$detection_rate == max_by_group
        ] <- "yes"

        # 更新原表
        genes_table$keep[
          match(dup_rows2$input, genes_table$input)
        ] <- dup_rows2$keep

        if (verbose) {
          n_resolved <- sum(dup_rows2$keep == "yes")
          cat("  Resolved", n_resolved, "duplicates by detection rate\n")
        }
      }
    }

    # 更新 conv_result
    conv_result <- genes_table

  } else {
    # Ortholog mode - simplified keep/duplication logic
    genes_table <- ortho_result
    genes_table$keep <- ifelse(genes_table$match_type == "found", "yes", "no")
    genes_table$duplication <- "no"
    ortho_result <- genes_table
  }

  # ==================== UPDATE FEATURES IF REQUESTED ====================

  keep_yes_idx <- genes_table$keep == "yes"
  n_keep_yes <- sum(keep_yes_idx, na.rm = TRUE)

  if (verbose) {
    cat("\n--- Keep Statistics ---\n")
    cat("  keep = yes:", n_keep_yes, "\n")
    cat("  keep = no:", sum(genes_table$keep == "no", na.rm = TRUE), "\n")
    cat("  keep = others:", sum(genes_table$keep == "others", na.rm = TRUE), "\n")
    cat("  duplication = yes:", sum(genes_table$duplication == "yes", na.rm = TRUE), "\n")
  }

  if (update_features) {
    if (verbose) cat("\nUpdating Seurat object features...\n")

    # 只保留 keep = "yes" 的基因
    if (n_keep_yes > 0) {
      # Get final symbols for keep = "yes" genes
      keep_final_symbols <- genes_table$final_symbol[keep_yes_idx]

      # Filter genes_table to only keep = "yes"
      if (sync_type == "convert") {
        conv_result <- genes_table[keep_yes_idx, , drop = FALSE]
      } else {
        ortho_result <- genes_table[keep_yes_idx, , drop = FALSE]
      }

      # Update new_features
      new_features <- keep_final_symbols

      # Get genes to keep (original names that map to keep = "yes")
      genes_to_keep <- orig_features[keep_yes_idx]

      # Get the default assay
      default_assay <- DefaultAssay(object)

      # Detect Seurat version by checking assay structure
      # V5 has @layers, V4 has @counts/@data/@scale.data
      assay_obj <- object[[default_assay]]
      seurat_v5 <- tryCatch({
        # V5 有 @layers 且非空，V4 有 @slots
        has_layers <- !is.null(assay_obj@layers) && length(assay_obj@layers) > 0
        # 额外检查：V5 的 @layers 是一个列表
        if (has_layers && is.list(assay_obj@layers)) {
          TRUE
        } else {
          FALSE
        }
      }, error = function(e) {
        FALSE
      })

      if (verbose) {
        if (seurat_v5) {
          cat("  Detected Seurat V5 object (has layers)\n")
        } else {
          cat("  Detected Seurat V4 object (has counts/data/scale.data slots)\n")
        }
      }

      if (seurat_v5) {
        # Seurat V5: use layer-based approach
        if (verbose) cat("  Processing Seurat V5 object...\n")

        # Get layer names
        layer_names <- names(assay_obj@layers)
        if (verbose) cat("  Found layers:", paste(layer_names, collapse = ", "), "\n")

        # 在 Seurat V5 中，需要确保 genes_to_keep 在 layer 的行名中存在
        # 先获取第一个 layer 的行名作为参考
        first_layer_name <- layer_names[1]
        first_layer <- assay_obj@layers[[first_layer_name]]
        layer_features <- rownames(first_layer)

        # 检查哪些 genes_to_keep 实际存在于 layer 中
        genes_to_keep_valid <- intersect(genes_to_keep, layer_features)
        if (length(genes_to_keep_valid) < length(genes_to_keep)) {
          if (verbose) {
            cat("  Note: Some genes not found in layers. Found:", length(genes_to_keep_valid),
                "of", length(genes_to_keep), "\n")
          }
        }

        # 调整 new_features 以匹配有效的基因
        valid_idx <- match(genes_to_keep_valid, genes_to_keep)
        new_features_valid <- new_features[valid_idx]

        # Create new layers list
        new_layers <- list()

        for (layer_name in layer_names) {
          if (verbose) cat("  Processing layer:", layer_name, "...\n")

          layer_data <- assay_obj@layers[[layer_name]]

          if (!is.null(layer_data)) {
            if (verbose) cat("    Layer class:", class(layer_data)[1], "\n")

            # Filter rows by genes_to_keep_valid
            tryCatch({
              # 检查 subset 后是否会有结果
              layer_features_in_layer <- intersect(genes_to_keep_valid, rownames(layer_data))

              if (length(layer_features_in_layer) == 0) {
                if (verbose) cat("    Warning: No matching genes found in this layer, skipping\n")
                next
              }

              if (verbose) cat("    Subsetting layer...\n")
              filtered_data <- layer_data[layer_features_in_layer, , drop = FALSE]

              # 检查 filtered_data 是否有效
              if (nrow(filtered_data) == 0 || ncol(filtered_data) == 0) {
                if (verbose) cat("    Warning: Subset resulted in empty matrix, skipping\n")
                next
              }

              # 调整 new_features_valid 以匹配实际 subset 的基因
              actual_idx <- match(layer_features_in_layer, genes_to_keep_valid)
              new_features_for_layer <- new_features_valid[actual_idx]

              if (verbose) cat("    Renaming rows...\n")
              rownames(filtered_data) <- new_features_for_layer

              new_layers[[layer_name]] <- filtered_data
              if (verbose) cat("    Layer done.\n")
            }, error = function(e) {
              if (verbose) cat("    Error in layer:", conditionMessage(e), "\n")
              # If subsetting fails, keep original
              new_layers[[layer_name]] <<- layer_data
            })
          }
        }

        # Update the assay
        if (verbose) cat("  Updating layers slot...\n")
        assay_obj@layers <- new_layers

        if (verbose) cat("  Updating features slot...\n")
        # **v3.2.2 修复**: 使用 LogMap() 构造 features，而不是 factor
        # 获取实际保留的基因名
        actual_kept_features <- unique(unlist(lapply(new_layers, rownames)))
        if (length(actual_kept_features) > 0) {
          assay_obj@features <- SeuratObject::LogMap(actual_kept_features)
        } else {
          warning("No features kept after subsetting. Features slot not updated.")
        }

        # Assign back to object
        object[[default_assay]] <- assay_obj

      } else {
        # Seurat V4: use subset then rename
        if (verbose) cat("  Processing Seurat V4 object...\n")

        # First, check which genes actually exist in the object
        current_features <- rownames(object)
        genes_found_idx <- genes_to_keep %in% current_features

        if (verbose) cat("  Genes to keep:", length(genes_to_keep), "\n")
        if (verbose) cat("  Genes found in object:", sum(genes_found_idx), "\n")

        # Subset the object
        object <- subset(x = object, features = genes_to_keep)

        # Get the actual features after subsetting
        actual_features <- rownames(object)
        n_actual <- length(actual_features)

        if (verbose) cat("  Features after subset:", n_actual, "\n")

        # Adjust new_features to match actual features
        # Find which of our genes_to_keep are in the actual object
        actual_idx <- match(actual_features, genes_to_keep)
        adjusted_new_features <- new_features[actual_idx]

        if (verbose) cat("  Adjusted new features:", length(adjusted_new_features), "\n")

        # Update gene names in assay slots
        assay_obj <- object[[default_assay]]

        # **v3.2.2 修复**: 确保 meta.features 与 data slots 一致
        # 先更新 slots 的行名
        # Update counts slot if exists
        if (!is.null(assay_obj@counts)) {
          if (nrow(assay_obj@counts) == length(adjusted_new_features)) {
            rownames(assay_obj@counts) <- adjusted_new_features
          }
        }
        # Update data slot if exists
        if (!is.null(assay_obj@data)) {
          if (nrow(assay_obj@data) == length(adjusted_new_features)) {
            rownames(assay_obj@data) <- adjusted_new_features
          }
        }
        # Update scale.data slot if exists
        if (!is.null(assay_obj@scale.data)) {
          if (nrow(assay_obj@scale.data) == length(adjusted_new_features)) {
            rownames(assay_obj@scale.data) <- adjusted_new_features
          }
        }

        # **v3.2.2 修复**: 更新 meta.features 以匹配 data slot 中的特征
        # 确保 meta.features 的行名与 data slot 一致
        if (!is.null(assay_obj@data)) {
          data_features <- rownames(assay_obj@data)
          meta_features <- assay_obj@meta.features
          # subset meta.features to match data_features
          common_features <- intersect(rownames(meta_features), data_features)
          if (length(common_features) > 0) {
            assay_obj@meta.features <- meta_features[common_features, , drop = FALSE]
          }
        }

        # Assign back
        object[[default_assay]] <- assay_obj
      }

      if (verbose) {
        cat("  Updated to", nrow(object), "features\n")
      }
    } else {
      warning("No genes with keep = 'yes' found. Object not updated.")
    }
  } else {
    if (verbose) cat("\nSkipping feature update (update_features = FALSE)\n")
    # Store the full genes_table for later use
    if (sync_type == "convert") {
      conv_result <- genes_table
    } else {
      ortho_result <- genes_table
    }
  }

  # ==================== STORE SYNC INFO ====================

  if (sync_type == "convert") {
    # Build sync_table with all columns
    sync_table <- data.frame(
      orig.feature = conv_result$input,
      new.feature = conv_result$final_symbol,
      symbol = conv_result$symbol,
      symbol_authority = conv_result$symbol_authority,
      gene_id = conv_result$gene_id,
      ensembl_id = conv_result$ensembl_id,
      chromosome = conv_result$chromosome,
      match_type = conv_result$match_type,
      final_symbol = conv_result$final_symbol,
      keep = conv_result$keep,
      duplication = conv_result$duplication,
      stringsAsFactors = FALSE
    )

    # Add n_final column if present
    if ("n_final" %in% names(conv_result)) {
      sync_table$n_final <- conv_result$n_final
    }

    summary_stats <- list(
      total_input = n_features,
      final_features = nrow(object),
      matched = sum(conv_result$match_type != "not_found", na.rm = TRUE),
      not_found = sum(conv_result$match_type == "not_found", na.rm = TRUE),
      keep_yes = sum(conv_result$keep == "yes", na.rm = TRUE),
      keep_no = sum(conv_result$keep == "no", na.rm = TRUE),
      keep_others = sum(conv_result$keep == "others", na.rm = TRUE),
      duplication_yes = sum(conv_result$duplication == "yes", na.rm = TRUE),
      update_features = update_features,
      match_type_counts = table(conv_result$match_type)
    )

    object@misc$gene_sync <- list(
      type = "convert",
      table = sync_table,
      summary = summary_stats,
      settings = list(
        species = species,
        from = from,
        to = to,
        symbol_type = symbol_type,
        update_features = update_features,
        timestamp = Sys.time()
      )
    )
  } else {
    # Ortholog sync
    sync_table <- data.frame(
      orig.feature = ortho_result$input,
      new.feature = ortho_result[[paste0(ortholog_to, "_symbol")]],
      keep = ortho_result$keep,
      duplication = ortho_result$duplication,
      stringsAsFactors = FALSE
    )

    # Add species-specific columns
    from_sym_col <- paste0(species, "_symbol")
    from_auth_col <- paste0(species, "_symbol_authority")
    from_id_col <- paste0(species, "_gene_id")
    to_sym_col <- paste0(ortholog_to, "_symbol")
    to_auth_col <- paste0(ortholog_to, "_symbol_authority")
    to_id_col <- paste0(ortholog_to, "_gene_id")

    sync_table$symbol <- ortho_result[[from_sym_col]]
    sync_table$symbol_authority <- ortho_result[[from_auth_col]]
    sync_table$gene_id <- ortho_result[[from_id_col]]
    sync_table[[paste0(ortholog_to, "_symbol")]] <- ortho_result[[to_sym_col]]
    sync_table[[paste0(ortholog_to, "_symbol_authority")]] <- ortho_result[[to_auth_col]]
    sync_table[[paste0(ortholog_to, "_gene_id")]] <- ortho_result[[to_id_col]]
    sync_table$match_type <- ortho_result$match_type

    summary_stats <- list(
      total_input = n_features,
      final_features = nrow(object),
      matched = sum(ortho_result$match_type == "found", na.rm = TRUE),
      not_found = sum(ortho_result$match_type == "not_found", na.rm = TRUE),
      keep_yes = sum(ortho_result$keep == "yes", na.rm = TRUE),
      keep_no = sum(ortho_result$keep == "no", na.rm = TRUE),
      keep_others = 0,
      duplication_yes = 0,
      update_features = update_features
    )

    object@misc$gene_sync <- list(
      type = "ortholog",
      table = sync_table,
      summary = summary_stats,
      settings = list(
        from_species = species,
        to_species = ortholog_to,
        symbol_type = symbol_type,
        update_features = update_features,
        timestamp = Sys.time()
      )
    )
  }

  if (verbose) {
    cat("\n--- Sync Summary ---\n")
    cat("Matched:", summary_stats$matched, "/", summary_stats$total_input, "\n")
    cat("Not found:", summary_stats$not_found, "/", summary_stats$total_input, "\n")
    cat("keep = yes:", summary_stats$keep_yes, "\n")
    cat("keep = no:", summary_stats$keep_no, "\n")
    if (summary_stats$keep_others > 0) {
      cat("keep = others:", summary_stats$keep_others, "\n")
    }
    cat("duplication = yes:", summary_stats$duplication_yes, "\n")
    cat("Final features:", summary_stats$final_features, "\n")
    cat("Time:", format(elapsed), "\n")
    cat("Done!\n")
  }

  return(object)
}

# 同时需要更新 get_gene_sync_table 和 print_gene_sync_summary 函数的文档
# 但函数本身不需要修改，因为它们已经可以正确处理 keep 和 duplication 列了
