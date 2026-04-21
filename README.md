# GSE146912 肾小球单细胞 —流程

本仓库提供 **统一配置 + Python 包 + 命令行入口**，用于从合并 `h5ad` 重跑预处理、聚类、注释、PEC 亚群与肾炎相关导出分析。
**转录因子调控用 SCENIC（推荐 pySCENIC）**，**细胞通讯用 CellChat（R）**；正式引用与 DOI 见 **[docs/METHODS_PAPER_ALIGNMENT.md](docs/METHODS_PAPER_ALIGNMENT.md)**。Python 中的 DoRothEA+ULM、LIANA 仅作探索性补充，默认不替代上述流程。

## 目录结构

```
config/                 # default.yaml、liu2023_kidney_int.yaml（论文对齐覆盖）、可选 local.yaml
docs/DATA_DOWNLOAD.md # 数据下载与输入文件说明
docs/METHODS_PAPER_ALIGNMENT.md  # 与 Kidney Int 2023 方法一致性与结果验收
src/gse146912_pipeline/  # 共用逻辑（注释、基因折叠、QC、PEC、通讯）
scripts/run_pipeline.py    # 主入口
results/                # 默认输出根（含 h5ad、figures、tables）
```

## 环境

```bash
cd /path/to/scRNA_Project
python -m venv .venv && source .venv/bin/activate
pip install -e ".[full]"    # scanpy + liana + decoupler + gseapy
pip install pybiomart scrublet   # 线粒体注释、双细胞
```

或使用 Conda（见 `environment.yml`，安装后仍需在仓库根目录执行 `pip install -e ".[full]"`）。

## 输入数据

将合并好的 **`GSE146912_Merged_Raw.h5ad`** 置于 `results/`（或修改 `config/default.yaml` 中 `paths.raw_h5ad`）。详见 **[docs/DATA_DOWNLOAD.md](docs/DATA_DOWNLOAD.md)**。

## 一键复现（推荐）

在**仓库根目录**执行（可通过环境变量 `SCRNA_PROJECT_ROOT` 指定根路径）：



```bash
# 仅主流程：QC（可配）→ 归一化判断 → HVG → 降维聚类 → marker 映射与打分 → 保存 analyzed h5ad
python scripts/run_pipeline.py --config config/default.yaml --stages main

# 默认配置 + liu2023 覆盖（仍须在 local.yaml 按原文补全 QC 等）
python scripts/run_pipeline.py --config config/default.yaml --local-config config/liu2023_kidney_int.yaml --stages main pec_subcluster

# 主流程 + PEC 亚群重聚类
python scripts/run_pipeline.py --config config/default.yaml --stages main pec_subcluster

# 肾炎对照子集：导出 symbol h5ad 供 R CellChat（config 中 injury_comm.enabled: true；默认不跑 LIANA）
python scripts/run_pipeline.py --config config/local.yaml --stages injury_comm

# 仅下游（pec_subcluster / injury_comm / 二者组合）：从上次 main 保存的 analyzed h5ad 读入，不重复 QC。
# 若缺少该文件会报错，请先运行 --stages main。
python scripts/run_pipeline.py --config config/default.yaml --stages pec_subcluster

# 全流程
python scripts/run_pipeline.py --config config/default.yaml --stages all
```

使用本地覆盖配置：

```bash
cp config/local.yaml.example config/local.yaml
python scripts/run_pipeline.py --config config/default.yaml --local-config config/local.yaml --stages main pec_subcluster
```

## 过程记录（provenance）

每次运行会在 **`results/GSE146912_run/provenance/`**（或你配置的 `output_run_dir`）写入：

- **`run_manifest.json`**：步骤列表、细胞数变化、关键参数、随机种子、marker 映射覆盖率摘要链接到分步文件。
- **`step_*.json`**：如 `batch_balance`、`marker_symbol_to_ensembl`、`injury_symbol_collapse` 等。

**基因 / marker**：`marker_symbol_to_ensembl` 中记录每个细胞类型的映射覆盖率与缺失 symbol；Ensembl→symbol 折叠见 `injury_symbol_collapse`（若运行通讯阶段）。

## 质量控制说明

在 `config/default.yaml` 的 `qc` 段配置：

- **线粒体比例**：`annotate_mito: true` 时通过 `scanpy.queries.mitochondrial_genes`（依赖 **pybiomart** 与网络）；失败时跳过线粒体比例，并在 manifest 中记录。
- **过滤**：`filter_min_genes` / `filter_max_mito_pct` 等，`null` 表示不启用。
- **双细胞**：`scrublet.enabled`；若数据已判定为 log 归一化（`xmax_skip_normalize` 启发式），将**自动跳过** Scrublet 并在 manifest 注明（log 矩阵上 Scrublet 不可靠）。
- **批次**：`batch_balance` 写出各 `batch` 细胞数及占比，过小将告警写入 JSON。

## 配置项摘要

| 配置路径 | 含义 |
|----------|------|
| `project.random_seed` | 全局随机种子 |
| `paths.raw_h5ad` / `output_run_dir` | 输入与输出目录 |
| `processing.batch_key` | HVG batch 键 |
| `processing.hvg_n_top_genes` | 高变基因数 |
| `leiden.resolution` | 全数据 Leiden 分辨率 |
| `annotation.marker_sets` | 细胞类型 marker（symbol） |
| `pec_subcluster.*` | PEC 内 HVG / Leiden 分辨率等 |
| `injury_comm.liana.enabled` | 默认 `false`；与论文 CellChat 主路径区分 |
| `injury_comm.liana.*` | 仅当 `enabled: true`：资源名、`expr_prop` |
| `qc.*` | 线粒体、过滤、Scrublet |



## 模块说明（`src/gse146912_pipeline`）

| 模块 | 职责 |
|------|------|
| `config.py` | 加载 YAML、合并 local、解析路径 |
| `provenance.py` | 输入指纹、manifest、分步 JSON |
| `markers.py` | marker symbol → Ensembl 与覆盖率 |
| `genes.py` | Ensembl → symbol 折叠（含 raw/X 维度不一致） |
| `injury.py` | Sample → injury_group |
| `qc.py` | 线粒体注释、QC 指标、过滤、批次统计、Scrublet |
| `annotation.py` | `score_genes` 与 `cell_type_major` |
| `pec.py` | PEC 内重聚类 |
| `injury_comm.py` | 对照/肾炎子集、symbol 折叠、导出 CellChat 用 h5ad；可选 LIANA |
| `pipeline_main.py` | 串联各步 |

## 引用

- **原始论文**：刘文斌（Liu WB）等，*Single cell landscape of parietal epithelial cells in healthy and diseased states*，*Kidney International* 2023；**GEO GSE146912**；doi:[10.1016/j.kint.2023.03.036](https://doi.org/10.1016/j.kint.2023.03.036)。
- 使用 **pySCENIC、CellChat、Scanpy、LIANA、decoupler** 等时，请分别按各工具官方文献引用。

详见 [docs/METHODS_PAPER_ALIGNMENT.md](docs/METHODS_PAPER_ALIGNMENT.md)。
