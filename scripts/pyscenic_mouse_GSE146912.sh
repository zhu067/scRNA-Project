#!/usr/bin/env bash
# pySCENIC（小鼠）：与刘文斌等 Kidney International 2023（GSE146912）原文 SCENIC 思路一致的主路径（regulon + AUCell）。
# DoRothEA/ULM 不能替代本流程作为主结论。
#
# 数据库从 https://resources.aertslab.org/cistarget/databases/mus_musculus/ 下载，
# 文件名随版本变化，请以页面为准（常见为 *.genes_vs_motifs.rankings.feather 与对应 annotations .tbl）。
#
# 依赖：按 pySCENIC 官方文档安装（示例）：
#   pip install pyscenic arboreto
#
# 必填环境变量：
#   SCENIC_EXPR      表达矩阵（tsv/csv/loom，行为细胞、列为基因 symbol）
#   SCENIC_OUT_DIR   输出目录
#   SCENIC_RANKING_FEATHER  cisTarget 排名库（至少一个 .feather，可多选追加参数）
#   SCENIC_ANNOT_TBL   motif 注释表（pyscenic ctx 的 --annotations_fname）
#
set -euo pipefail

: "${SCENIC_EXPR:?}"
: "${SCENIC_OUT_DIR:?}"
: "${SCENIC_RANKING_FEATHER:?}"
: "${SCENIC_ANNOT_TBL:?}"

mkdir -p "$SCENIC_OUT_DIR"
ADJ="$SCENIC_OUT_DIR/adjacencies.tsv"
REG="$SCENIC_OUT_DIR/regulons.csv"
AUC="$SCENIC_OUT_DIR/auc_mtx.csv"

echo "== (1/3) GRNBoost2：邻接矩阵 =="
pyscenic grn "$SCENIC_EXPR" "$ADJ" --num_workers 8

echo "== (2/3) cisTarget：剪枝 regulons =="
# 若有多个 feather 库，可在此行追加更多路径
pyscenic ctx "$ADJ" "$SCENIC_RANKING_FEATHER" --annotations_fname "$SCENIC_ANNOT_TBL" -o "$REG" --num_workers 8

echo "== (3/3) AUCell：regulon 活性 =="
pyscenic aucell "$SCENIC_EXPR" "$REG" -o "$AUC"

echo "完成。请记录：pyscenic 版本、上述四个输入路径、数据库文件名，并写入 provenance。"
