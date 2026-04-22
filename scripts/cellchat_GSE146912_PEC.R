# CellChat（R）：与刘文斌等 Kidney International 2023（GSE146912, PMID 37100348；
# 「肾小球壁层上皮细胞亚群与新月体形成机制及治疗」相关单细胞研究）中
# glomerular cell–cell communication 分析方向一致的主路径。
# Python 管线中的 LIANA 为可选补充，不等价于原文 CellChat 结果。
#
# 自动对照 vs 肾炎：按 meta$injury_group 拆成两个 CellChat，merge 后出差异图与矩阵。
#
# 依赖示例：
#   BiocManager::install(c("zellkonverter", "SingleCellExperiment", "ggplot2"))
#   devtools::install_github("sqjin/CellChat")
#
# 请先运行 Python：python scripts/run_pipeline.py --stages main injury_comm
# 并确认存在 Python 导出的计数表（供样本量检查）：
#   results/GSE146912_run/tables/PEC_injury_comm/injury_group_cell_type_major_counts.csv
#
# 在项目根目录： Rscript scripts/cellchat_GSE146912_PEC.R

suppressPackageStartupMessages({
  if (!requireNamespace("zellkonverter", quietly = TRUE))
    stop("请安装: BiocManager::install('zellkonverter')")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("请安装: BiocManager::install('SingleCellExperiment')")
  if (!requireNamespace("CellChat", quietly = TRUE))
    stop("请安装 CellChat: https://github.com/sqjin/CellChat")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("请安装 ggplot2（CellChat 作图常用）")
})

MIN_CELLS_PER_GROUP_TYPE <- 10L

h5 <- "results/GSE146912_run/GSE146912_for_CellChat_symbol_injury.h5ad"
counts_csv <- "results/GSE146912_run/tables/PEC_injury_comm/injury_group_cell_type_major_counts.csv"
out_base <- "results/GSE146912_run/tables/PEC_injury_comm/cellchat_R_diff"

if (!file.exists(h5)) stop("未找到 ", h5)

assert_injury_crosstab <- function(path, min_cells = MIN_CELLS_PER_GROUP_TYPE) {
  if (!file.exists(path)) {
    stop(
      "未找到计数表 ", path,
      "。请先跑 Python injury_comm 阶段以导出 injury_group x cell_type_major。"
    )
  }
  tab <- read.csv(path, row.names = 1, check.names = FALSE)
  rn <- rownames(tab)
  if (!all(c("control", "nephritis") %in% rn)) {
    stop("计数表行须包含 injury_group：control、nephritis（当前: ", paste(rn, collapse = ", "), "）")
  }
  bad <- character(0)
  for (ct in colnames(tab)) {
    c_ctrl <- suppressWarnings(as.integer(tab["control", ct]))
    c_ne <- suppressWarnings(as.integer(tab["nephritis", ct]))
    if (is.na(c_ctrl)) c_ctrl <- 0L
    if (is.na(c_ne)) c_ne <- 0L
    if (c_ctrl > 0L && c_ne > 0L && (c_ctrl < min_cells || c_ne < min_cells)) {
      bad <- c(bad, sprintf("%s (control=%d, nephritis=%d)", ct, c_ctrl, c_ne))
    }
  }
  if (length(bad)) {
    stop(
      "以下 cell_type_major 在两组均存在但样本量 < ", min_cells, "：\n",
      paste(bad, collapse = "\n"),
      "\n请放宽过滤或提高测序深度后再跑 CellChat。"
    )
  }
  invisible(tab)
}

#' 自 createCellChat 起跑完整 CellChat 流程至 aggregateNet（单 injury 组）
run_one_group <- function(sce, group_value, celltype_levels) {
  meta <- as.data.frame(SingleCellExperiment::colData(sce))
  keep <- !is.na(meta$injury_group) & meta$injury_group == group_value
  if (!any(keep)) stop("injury_group==", group_value, " 的细胞数为 0")
  sce_g <- sce[, keep]
  data.input <- SummarizedExperiment::assay(sce_g)
  meta_g <- as.data.frame(SingleCellExperiment::colData(sce_g))
  meta_g$cell_type_major <- factor(as.character(meta_g$cell_type_major), levels = celltype_levels)

  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta_g,
    group.by = "cell_type_major"
  )
  cellchat@DB <- CellChat::CellChatDB.mouse
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", trim = 0.1)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat
}

write_net_csvs <- function(cellchat, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cellchat@net$count, file.path(outdir, "cellchat_interaction_count_matrix.csv"))
  utils::write.csv(cellchat@net$weight, file.path(outdir, "cellchat_interaction_weight_matrix.csv"))
}

netp_prob_matrix <- function(cc) {
  np <- cc@netP
  if (is.null(np)) return(NULL)
  p <- if (is.list(np)) np[["prob"]] else tryCatch(methods::slot(np, "prob"), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  if (inherits(p, "dgCMatrix")) p <- as.matrix(p)
  if (!is.matrix(p)) return(NULL)
  p
}

safe_pathway_diff_csv <- function(cc_ctrl, cc_ne, path_out) {
  ok <- FALSE
  tryCatch(
    {
      p1 <- netp_prob_matrix(cc_ctrl)
      p2 <- netp_prob_matrix(cc_ne)
      if (is.null(p1) || is.null(p2)) return(invisible(FALSE))
      if (!all(dim(p1) == dim(p2))) return(invisible(FALSE))
      dp <- p2 - p1
      rownames(dp) <- rownames(p1)
      colnames(dp) <- colnames(p1)
      utils::write.csv(dp, path_out)
      ok <- TRUE
    },
    error = function(e) {
      message("pathway 差异矩阵跳过: ", conditionMessage(e))
    }
  )
  invisible(ok)
}

# --- main ---
assert_injury_crosstab(counts_csv, MIN_CELLS_PER_GROUP_TYPE)

sce <- zellkonverter::readH5AD(h5)
meta_full <- as.data.frame(SingleCellExperiment::colData(sce))
if (!"injury_group" %in% colnames(meta_full)) stop("colData 缺少 injury_group")
if (!"cell_type_major" %in% colnames(meta_full)) stop("colData 缺少 cell_type_major")

celltype_levels <- sort(unique(as.character(meta_full$cell_type_major)))

cc_control <- run_one_group(sce, "control", celltype_levels)
cc_nephritis <- run_one_group(sce, "nephritis", celltype_levels)

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
write_net_csvs(cc_control, file.path(out_base, "control"))
write_net_csvs(cc_nephritis, file.path(out_base, "nephritis"))

diff_dir <- file.path(out_base, "diff")
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)

diff_count <- cc_nephritis@net$count - cc_control@net$count
diff_weight <- cc_nephritis@net$weight - cc_control@net$weight
utils::write.csv(diff_count, file.path(diff_dir, "cellchat_diff_count_matrix.csv"))
utils::write.csv(diff_weight, file.path(diff_dir, "cellchat_diff_weight_matrix.csv"))

pathway_out <- file.path(diff_dir, "cellchat_diff_pathway_prob.csv")
had_pathway <- safe_pathway_diff_csv(cc_control, cc_nephritis, pathway_out)
if (!had_pathway) {
  message("未写出 pathway 概率差异（当前 CellChat 对象 netP@prob 不可用或维度不一致）。")
}

object.list <- list(control = cc_control, nephritis = cc_nephritis)
cellchat_merged <- CellChat::mergeCellChat(object.list, add.names = names(object.list))

grDevices::pdf(file.path(diff_dir, "cellchat_diff_interaction_count.pdf"), width = 9, height = 7)
CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count")
grDevices::dev.off()

grDevices::pdf(file.path(diff_dir, "cellchat_diff_interaction_weight.pdf"), width = 9, height = 7)
CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
grDevices::dev.off()

if (requireNamespace("ggplot2", quietly = TRUE)) {
  tryCatch(
    {
      p1 <- CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count")
      ggplot2::ggsave(
        file.path(diff_dir, "cellchat_diff_interaction_count.png"),
        p1,
        width = 9,
        height = 7,
        dpi = 150
      )
    },
    error = function(e) message("diffInteraction count 导出 PNG 跳过: ", conditionMessage(e))
  )
  tryCatch(
    {
      p2 <- CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
      ggplot2::ggsave(
        file.path(diff_dir, "cellchat_diff_interaction_weight.png"),
        p2,
        width = 9,
        height = 7,
        dpi = 150
      )
    },
    error = function(e) message("diffInteraction weight 导出 PNG 跳过: ", conditionMessage(e))
  )
}

message("完成。输出目录: ", normalizePath(out_base))
message("差异矩阵（nephritis − control）与 diffInteraction 图位于: ", normalizePath(diff_dir))
