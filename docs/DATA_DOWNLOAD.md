# 数据获取说明（GSE146912）

##  accession

- **GEO**: [GSE146912](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146912)
- 数据：小鼠肾小球单细胞，含健康与多种损伤模型样本。

## 本仓库期望的输入文件

默认配置（`config/default.yaml`）假定以下文件已就绪：

| 路径（相对仓库根） | 说明 |
|-------------------|------|
| `results/GSE146912_Merged_Raw.h5ad` | 合并后的 AnnData；`obs` 至少含 `Sample`、`batch`；`var` 为小鼠 **Ensembl** 基因 ID |

请将预处理得到的合并矩阵放到该路径，或修改配置中的 `paths.raw_h5ad`。

## 从 GEO 自行下载时的提示

1. 在 GEO 页面下载 **Series Matrix** 或作者提供的 **processed** 数据；原始 `mtx`/`barcodes`/`genes` 需自行用 Scanpy/Seurat 读入并合并。
2. 若仅下载 **10X** 格式 per-sample，需循环读取并 `concatenate`，统一基因（Ensembl 或 symbol）后再写出 `h5ad`。
3. 合并后建议记录：每样本细胞数、基因数、是否与本文使用的 `batch` 编码一致。

## 公开镜像（可选）

- Broad Single Cell Portal 等可能提供同一数据集的处理版本；使用第三方文件时请在 `provenance/run_manifest.json` 中核对 `input_h5ad` 指纹与样本元数据。

## 模型说明

- 样本名中含 `nephritis` 的模型在文献中为 **NTS（肾毒性血清肾炎）**，与临床 **抗 GBM** 不同；比较分组见配置 `injury_comm` 中的 token。

## 分析入口

- 方法学、结果验收与主流程（SCENIC + CellChat）：**[METHODS_PAPER_ALIGNMENT.md](METHODS_PAPER_ALIGNMENT.md)**。
- 命令行推荐覆盖层：**`config/liu2023_kidney_int.yaml`**（与 `default.yaml` 联用；QC 等仍按原文在 `local.yaml` 补全）。
