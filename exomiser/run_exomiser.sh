#!/bin/bash
# Exomiser的运行器
# pzw
# 20250728
# 一并封装到exomiser的镜像里，因此调用的是镜像中的路径

# 默认值
DEFAULT_YAML="/exomiser-cli-14.1.0/config/default.yaml"
RUNNING_DIR="/run_data"
APPLICATION_PROPERTIES="/exomiser-cli-14.1.0/application.properties"

# 显示帮助信息
show_help() {
    cat << EOF
用法: $0 [选项]

必需参数:
  --sample_id SAMPLE_ID     样本ID
  --vcf_file VCF_FILE       VCF文件路径
  --output_dir OUTPUT_DIR   输出目录

可选参数:
  --genome GENOME           基因组版本 (如: hg19, hg38)
  --hpo_string HPO_STRING   HPO术语，用逗号分隔
  --ped_file PED_FILE       PED文件路径
  --data_dir DATA_DIR       Exomiser数据目录路径
  --data_version VERSION    数据库版本号 (phenotype, hg19, hg38统一版本)
  -h, --help               显示此帮助信息

示例:
  $0 --sample_id sample1 --vcf_file input.vcf --output_dir /path/to/output --genome hg38 --hpo_string "HP:0001250,HP:0000252" --data_dir /data/exomiser --data_version 2409
EOF
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample_id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --vcf_file)
            VCF_FILE="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --hpo_string)
            HPO_STRING="$2"
            shift 2
            ;;
        --ped_file)
            PED_FILE="$2"
            shift 2
            ;;
        --data_dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --data_version)
            DATA_VERSION="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "未知参数: $1"
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$SAMPLE_ID" || -z "$VCF_FILE" || -z "$OUTPUT_DIR" ]]; then
    echo "错误: 缺少必需参数"
    show_help
    exit 1
fi

# 检查yq是否安装
if ! command -v yq &> /dev/null; then
    echo "错误: 需要安装yq工具来处理YAML文件"
    echo "请安装yq: https://github.com/mikefarah/yq"
    exit 1
fi

# 更新application.properties文件的函数
update_application_properties() {
    local prop_file="$1"
    local key="$2"
    local value="$3"
    
    # 检查属性是否存在（包括被注释的）
    if grep -q "^#*${key}=" "$prop_file"; then
        # 如果存在，更新值（取消注释并设置新值）
        sed -i "s|^#*${key}=.*|${key}=${value}|" "$prop_file"
    else
        # 如果不存在，添加新属性
        echo "${key}=${value}" >> "$prop_file"
    fi
}

# 创建运行目录
mkdir -p "$RUNNING_DIR"

# 获取绝对路径
VCF_FILE=$(realpath "$VCF_FILE")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

# 生成配置文件路径
CONFIG_YAML="$RUNNING_DIR/${SAMPLE_ID}_exomiser_config.yaml"

# 复制默认配置文件
cp "$DEFAULT_YAML" "$CONFIG_YAML"

# 更新VCF文件路径
yq eval ".analysis.vcf = \"$VCF_FILE\"" -i "$CONFIG_YAML"

# 更新基因组版本（如果提供）
if [[ -n "$GENOME" ]]; then
    yq eval ".analysis.genomeAssembly = \"$GENOME\"" -i "$CONFIG_YAML"
fi

# 更新HPO术语（如果提供）
if [[ -n "$HPO_STRING" ]]; then
    # 将逗号分隔的字符串转换为YAML数组
    IFS=',' read -ra HPO_ARRAY <<< "$HPO_STRING"
    yq eval ".analysis.hpoIds = []" -i "$CONFIG_YAML"
    for hpo in "${HPO_ARRAY[@]}"; do
        hpo=$(echo "$hpo" | xargs)  # 去除空格
        yq eval ".analysis.hpoIds += [\"$hpo\"]" -i "$CONFIG_YAML"
    done
fi

# 更新PED文件（如果提供）
if [[ -n "$PED_FILE" ]]; then
    PED_FILE=$(realpath "$PED_FILE")
    yq eval ".analysis.ped = \"$PED_FILE\"" -i "$CONFIG_YAML"
    yq eval ".analysis.proband = \"$SAMPLE_ID\"" -i "$CONFIG_YAML"
fi

echo "配置文件已生成: $CONFIG_YAML"

# 更新application.properties文件（如果提供了相关参数）
if [[ -n "$DATA_DIR" || -n "$DATA_VERSION" ]]; then
    echo "更新application.properties文件..."
    
    # 备份原始文件
    cp "$APPLICATION_PROPERTIES" "${APPLICATION_PROPERTIES}.backup"
    
    # 更新数据目录
    if [[ -n "$DATA_DIR" ]]; then
        update_application_properties "$APPLICATION_PROPERTIES" "exomiser.data-directory" "$DATA_DIR"
        echo "已设置数据目录: $DATA_DIR"
    fi
    
    # 更新数据版本
    if [[ -n "$DATA_VERSION" ]]; then
        update_application_properties "$APPLICATION_PROPERTIES" "exomiser.phenotype.data-version" "$DATA_VERSION"
        update_application_properties "$APPLICATION_PROPERTIES" "exomiser.hg19.data-version" "$DATA_VERSION"
        update_application_properties "$APPLICATION_PROPERTIES" "exomiser.hg38.data-version" "$DATA_VERSION"
        echo "已设置数据版本: $DATA_VERSION"
    fi
    
    echo "application.properties文件已更新"
fi

# 运行Exomiser
echo "开始运行Exomiser..."
java -jar /exomiser-cli-14.1.0/exomiser-cli-14.1.0.jar \
    --analysis "$CONFIG_YAML" \
    --output-directory "$OUTPUT_DIR" \
    --output-filename "$SAMPLE_ID" \
    --spring.config.location=/exomiser-cli-14.1.0/application.properties

echo "Exomiser运行完成"
