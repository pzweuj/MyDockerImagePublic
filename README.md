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

CNVkit 改版增加了参数调整脚本，用于优化 CNV 检测精度。

> ⚠️ **注意**：官方不建议调整该参数，但当进行高分辨率的 CNV 检测时，较低/较高的 GC 比例设定可能会导致一些真实 CNV 区域被过滤掉，请谨慎进行调整。

**调整 GC 比例参数示例：**
```bash
python /opt/conda/bin/cnvkit_params_modify.py --force_rewrite True --GC_MIN_FRACTION 0.25
```

**调整自动检测性别参数：**
将默认按 antitarget 检测修改为默认按 target 检测
```bash
python /opt/conda/bin/cnvkit_params_modify.py --reference_auto_model True
```

**Singularity 使用示例：**
使用 `exec` 运行时，需加入 `--writable-tmpfs` 参数
```bash
singularity exec --writable-tmpfs cnvkit_v0.9.11.p4.sif bash -c \
  "python /opt/conda/bin/cnvkit_params_modify.py --reference_auto_model True && \
   cnvkit.py reference coverage/*.{,anti}targetcoverage.cnn \
   --fasta human_g1k_v37_decoy.fasta -o reference.cnn"
```

