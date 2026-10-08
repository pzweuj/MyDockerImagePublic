# MyDockerImagePublic

用于部署 Docker 镜像到 GitHub Container Registry (GHCR)。

📦 [我的镜像库](https://github.com/pzweuj?tab=packages)

## 使用方法

### 1. 创建 GHCR Token

1. 点击 GitHub 头像 → **Settings**
2. 选择左下角的 **Developer settings**
3. 选择 **Personal access tokens** → **Tokens (classic)**
4. 点击 **Generate new token** → 选择 **classic**
5. 输入描述，选择 **repo** 权限，时间最长选择一年
6. 点击 **Generate token**，复制生成的 token（需要拥有读写权限）

### 2. 配置仓库密钥

1. 在 GitHub 仓库的 **Settings** 中，选择 **Secrets and variables** → **Actions**
2. 点击 **New repository secret**
3. Name: `GHCR_PAT`
4. Value: 刚刚复制的 token
5. 点击 **Add secret**

### 3. 部署镜像

1. 在仓库中为不同项目创建不同的文件夹
2. 放入 `Dockerfile`
3. 输入 commit message（此时 commit message 即为部署后的 tag）
4. 点击 **commit and push**，等待部署完成

## 生物信息镜像

| 镜像名称 | 功能描述 | 包含工具/特性 | 应用场景 |
|---------|----------|---------------|----------|
| **Mapping** | 序列比对和质控工具集 | fastp、bwa、samtools、sambamba、bamdst | 测序数据预处理和比对 |
| **Optitype** | HLA 分型检测软件 | HLA 分型算法 | 人类白细胞抗原分型分析 |
| **CNVkit** | 拷贝数变异检测工具 | CNVkit + 参数调整脚本 | 高分辨率 CNV 检测 |
| **AutoMap** | ROH 检测软件 | 同源性区段分析 | 全外显子测序数据分析 |
| **AutoCNV** | CNV 分析工具 | AutoCNV 封装版本 | 拷贝数变异分析（部分数据库缺失） |
| **Whatshap** | 单倍型计算软件 | WhatsHap、tabix、bgzip | 单倍型相位分析 |
| **Manta** | 结构变异检测工具 | Illumina SV 检测算法 | 大片段结构变异识别 |
| **Exomiser** | 有害性预测工具 | ACMG 标准评估 | 变异致病性评估 |
| **MSIsensor-pro** | 微卫星不稳定性分析 | MSI 检测算法 | NGS 数据 MSI 状态评估 |

### CNVkit 参数调整说明

镜像基座为 CNVkit v0.9.14。参数脚本只修改 `cnvlib/params.py` 里的六个常量，不覆盖文件其余内容，也不再修改 `reference.py`。

> 官方不建议调整这些常量。高分辨率检测时，GC 比例过窄可能滤掉真实 CNV 区域，请先用固定 BAM 和 BED 对比后再用于生产。

**查看当前生效值：**

```bash
python /opt/conda/bin/cnvkit_params_modify.py
```

**调整 GC 比例：**

```bash
python /opt/conda/bin/cnvkit_params_modify.py --force_rewrite True --GC_MIN_FRACTION 0.25
```

`--force_rewrite True` 可以省略。`--force_rewrite False` 不会写文件。

**性别推断：**

v0.9.14 在 target 与 antitarget 冲突时比较 chrX 证据强度，捕获 panel 上通常采用 target。`--reference_auto_model` 仍可传入，但只打印警告，不再改源码。

用旧补丁生成过的 Reference 不会随镜像升级自动更正。生产 panel 应使用 v0.9.14 重新构建 Reference，并用同一批 BAM、BED 做对比。差异若集中在 chrX/chrY，优先核对性别推断。

**Singularity：**

`exec` 写参数文件时需要 `--writable-tmpfs`。

```bash
singularity exec --writable-tmpfs cnvkit_v0.9.14.sif bash -c \
  "python /opt/conda/bin/cnvkit_params_modify.py --GC_MIN_FRACTION 0.25 && \
   cnvkit.py reference coverage/*.{,anti}targetcoverage.cnn \
   --fasta human_g1k_v37_decoy.fasta -o reference.cnn"
```

