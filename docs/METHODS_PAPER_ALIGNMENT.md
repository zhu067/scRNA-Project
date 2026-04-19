# 与原始论文方法对齐（GSE146912）

## 论文信息（您所指原文）

**课题常用中文表述**：肾小球壁层上皮细胞亚群组成及新月体形成的机制与治疗（**第一作者：刘文斌**；团队含 Xu AL 等）。

**期刊正式英文题名**：*Single cell landscape of parietal epithelial cells in healthy and diseased states*（Liu WB *et al.*）。

| 项目 | 内容 |
|------|------|
| 期刊 | *Kidney International* |
| 年份 / 卷期 | 2023; 104(1):108–123（以出版社页码为准） |
| DOI | [10.1016/j.kint.2023.03.036](https://doi.org/10.1016/j.kint.2023.03.036) |
| PMID | [37100348](https://pubmed.ncbi.nlm.nih.gov/37100348/) |
| 数据 | GEO [**GSE146912**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146912) |

下文「原文」「论文」均指上述 **刘文斌等 / Kidney International 2023** 这项研究。

---

## 结果优先：与论文一致的最小路径

若要以**可发表的、与原文可比的结果**为目标（而非仅「跑通」），请按下面顺序自检；**任何替代流程只能作附录或敏感性分析，不能顶替主结论**。

1. **输入矩阵**  
   - 使用与 GEO GSE146912 及原文预处理一致的 **`h5ad`**（基因 ID、是否已合并批次、是否为**原始计数**）。  
   - 本仓库主流程若根据 `X` 最大值判断为「已 log/缩放」而**跳过** `normalize_total` + `log1p`，会在 `provenance/` 写入 **`paper_replication_warning`**：此时与原文 Scanpy 归一化链可能不一致，且 **SCENIC/GRNBoost 通常需要按 pySCENIC 文档准备表达矩阵**，请勿混用未声明的变换。

2. **基础分析（Scanpy）**  
   - 论文主文基于 **Scanpy** 类流程；具体 **QC 阈值、HVG 数、Leiden 分辨率、批次校正与否** 必须以 **原文及 Supplementary** 为准，在 `config/liu2023_kidney_int.yaml` 或 `config/local.yaml` 中**显式写出**，并在 `run_manifest.json` 中可追溯。

3. **主图级下游**  
   - **TF / regulon：SCENIC（或 pySCENIC）**，数据库版本与参数记入 provenance。  
   - **细胞通讯：CellChat（R）**，小鼠库版本与 `injury_group` 分层策略与原文一致；Python **LIANA 默认关闭**。

4. **结果验收（建议）**  
   - 与原文或补充表对照：**总细胞数、主要类型占比、PEC 亚群数量与命名逻辑、关键 marker**；通讯与 TF 结果只做定性方向一致时，需检查输入与参数是否与 Methods 一致。  
   - 若无法复现，优先排查：**输入是否同源、归一化是否被跳过、基因 ID 是否一致**，其次才是随机种子与软件小版本。

---

## 1. 转录因子活性：SCENIC（非 DoRothEA+ULM）

- 原文使用 **SCENIC**（single-cell regulatory network inference and clustering；Aibar *et al.*, *Nat Methods* 2017）推断 **regulon** 并在细胞中评估活性。
- **decoupler + DoRothEA + ULM** 属于基于先验 TF–靶基因网络的统计推断，**算法与假设与 SCENIC 不同**，不能作为与论文 Figure 一一对应的复现路径。
- 推荐在独立环境中运行 **pySCENIC**（或 R SCENIC），步骤与数据库版本需自行与原文 / 补充材料核对。本仓库提供命令模板：**`scripts/pyscenic_mouse_GSE146912.sh`**（需下载 cisTarget 排名数据库与 motif 库，见脚本内注释与 [aertslab 资源页](https://resources.aertslab.org/cistarget/)）。

### 输入建议

- 使用与本流程一致的 **`GSE146912_merged_analyzed.h5ad`** 中 **PEC 子集**或全文，导出为 pySCENIC 接受的表达矩阵格式（loom / csv），基因名为 **symbol** 或与所用数据库一致的 ID。
- 记录 **GRNBoost2 / co-expression**、**cisTarget** 数据库文件名与版本、**AUC** 阈值，并写入 `provenance/`（与主流程 manifest 并列说明即可）。

## 2. 细胞通讯：CellChat（R）

- 原文对信号网络的解读基于 **CellChat**（Jin *et al.*, *Nat Commun* 2021）等 **cell–cell communication** 分析框架（见论文 Methods 中 dynamic signaling / CellChat 相关描述）。
- 本仓库在 `injury_comm` 阶段默认 **仅导出** `GSE146912_for_CellChat_symbol_injury.h5ad`，主分析路径为 R 脚本 **`scripts/cellchat_GSE146912_PEC.R`**（小鼠数据库 `CellChatDB.mouse`，流程与论文方向一致；组间比较需按 `injury_group` 分层构建两个 CellChat 对象后做差异可视化）。
- **LIANA**（Python）整合多数据库与多种打分，可作为 **补充探索**；配置项 `injury_comm.liana.enabled` 默认为 `false`，避免与「论文主流程」混淆。

## 3. 本仓库中「探索性」分析的定位

| 组件 | 论文主流程 | 本仓库默认 |
|------|------------|------------|
| TF / regulon | SCENIC | 文档 + `pyscenic_mouse_GSE146912.sh`，**不**默认跑 DoRothEA |
| 细胞通讯 | CellChat（R） | 导出 h5ad + `cellchat_GSE146912_PEC.R` |
| LIANA | 非原文必需 | `liana.enabled: false` |
| DoRothEA + ULM | 非原文方法 | 笔记本中可选，`EXPLORATORY_PYTHON_MULTIMETHOD = True` 时运行 |

## 4. 其它使用同一数据集的文献

若你复现的是后续引用 GSE146912 的文章（例如部分工作流写明 **pySCENIC 0.12.x** + **CellChat 1.6.x**），请以**该文 Methods** 为准，并将软件版本与数据库文件记入 `provenance`。

## 5. 引用

- **原始研究**：Liu WB *et al.* *Single cell landscape of parietal epithelial cells in healthy and diseased states.* *Kidney International* 2023;104(1):108–123. doi:10.1016/j.kint.2023.03.036（GEO GSE146912；中文主题见上文）。
- SCENIC：Aibar *et al.*, *Nat Methods* 2017。
- pySCENIC：Van de Sande *et al.*, *Nat Protoc* 2020（及所用版本说明）。
- CellChat：Jin *et al.*, *Nat Commun* 2021。
