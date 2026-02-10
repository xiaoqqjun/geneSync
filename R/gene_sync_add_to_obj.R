# gene_sync_add_to_obj 函数更新版本
# 主要更改：
# 1. 添加 update_features 参数（默认 TRUE）
# 2. 添加 keep 列（yes/no/others）
# 3. 添加 duplication 列（yes/no）
# 4. 根据检测率处理重复基因

# 请将下面的代码复制并替换 gene_sync_add_to_obj.R 中的整个函数

gene_sync_add_to_obj <- function(object,
                                  species = c("homo", "mmu", "rat"),
                                  symbol_type = c("final_symbol", "symbol_authority", "symbol"),
                                  ortholog_to = NULL,
                                  update_features = TRUE,
                                  verbose = TRUE) {

  species <- match.arg(species)
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
    conv_result <- gene_convert(orig_features, species = species)
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
        # 使用 GetAssayData 获取表达数据
        expr_mat <- Seurat::GetAssayData(object, slot = "data")[dup_inputs, , drop = FALSE]
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

      # Filter the Seurat object
      object <- object[keep_yes_idx, , drop = FALSE]

      # Update rownames
      seurat_v5 <- tryCatch(
        packageVersion("Seurat") >= "5.0.0",
        error = function(e) FALSE
      )

      if (seurat_v5) {
        if ("RNA" %in% names(object@assays)) {
          object@assays$RNA@features <- new_features
        }
        rownames(object) <- new_features
      } else {
        rownames(object) <- new_features
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
