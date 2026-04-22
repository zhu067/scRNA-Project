# CellChat（R）
# 「肾小球壁层上皮细胞亚群与新月体形成机制及治疗」相关单细胞研究）中
# Python 管线中的 LIANA 为可选补充。
#
# 对照 vs 肾炎：分别对 meta$injury_group == "control" / "nephritis" 子集构建 CellChat，
# 再使用 CellChat::netVisual_diffInteraction 等函数比较（见文末提示）。
#
# 依赖示例：
#   BiocManager::install(c("zellkonverter", "SingleCellExperiment"))
#   devtools::install_github("sqjin/CellChat")   # 或按 CellChat 官方说明安装
#
# 请先运行 Python 笔记本导出：
#   results/GSE146912_run/GSE146912_for_CellChat_symbol_injury.h5ad
#
# 在项目根目录： source("scripts/cellchat_GSE146912_PEC.R")

suppressPackageStartupMessages({
  if (!requireNamespace("zellkonverter", quietly = TRUE))
    stop("请安装: BiocManager::install('zellkonverter')")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("请安装: BiocManager::install('SingleCellExperiment')")
  if (!requireNamespace("CellChat", quietly = TRUE))
    stop("请安装 CellChat: https://github.com/sqjin/CellChat")
})

h5 <- "results/GSE146912_run/GSE146912_for_CellChat_symbol_injury.h5ad"
if (!file.exists(h5)) stop("未找到 ", h5)

sce <- zellkonverter::readH5AD(h5)
# SCE：行=基因，列=细胞（与 zellkonverter 对 h5ad 的约定一致）
data.input <- SummarizedExperiment::assay(sce)
meta <- as.data.frame(SingleCellExperiment::colData(sce))

cellchat <- CellChat::createCellChat(
  object = data.input,
  meta = meta,
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

outdir <- "results/GSE146912_run/tables/PEC_injury_comm/cellchat_R"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(cellchat@net$count, file.path(outdir, "cellchat_interaction_count_matrix.csv"))
write.csv(cellchat@net$weight, file.path(outdir, "cellchat_interaction_weight_matrix.csv"))

message("完成。矩阵已写入 ", normalizePath(outdir))
message("差异比较：对 meta$injury_group 分别子集重复本流程，使用 CellChat::netVisual_diffInteraction 对比两个 cellchat 对象。")
