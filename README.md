# MyDockerImagePublic

用于部署docker镜像到ghcr。

[我的镜像库](https://github.com/pzweuj?tab=packages)

## 使用方法

1，创建一个自己的Ghcr Token，点击自己的头像，选择Setting，选择左下角的Developer settings，选择Personal access tokens，点击Tokens(classic)，点击Generate new token,同样选择classic，输入描述，选择repo，时间最长是选择一年，点击Generate token，注意token需要拥有读写权限，复制token；

2，在github的setting中，选择Secrets，点击New repository secret，输入Name为GHCR_PAT，输入Value为刚刚复制的token，点击Add secret；

3，在仓库中，对不同项目新建不同的文件夹，放入dockerfile，输入commit message，此时commit message即是部署后的tag，点击commit and push，等待部署完成。

## 生物信息镜像

### Mapping

包含fastp、bwa、samtools、sambamba、bamdst。

### Optitype

HLA分型检测软件。

### CNVkit

cnvkit改版，增加一个脚本用于调整params参数。注意，官方不建议调整该参数，但当进行高分辨率的CNV检测时，较低/较高的GC比例设定可能会导致一些真实CNV区域被过滤掉，谨慎进行调整。

示例
```bash
python /opt/conda/bin/cnvkit_params_modify.py --force_rewrite True --GC_MIN_FRACTION 0.25
```

### AutoMap

用于WES的ROH检测软件。

### AutoCNV

AutoCNV的封装版本，部分数据库缺失。

### Whatshap

单倍型计算软件WhatsHap镜像，包含tabix、bgzip。

