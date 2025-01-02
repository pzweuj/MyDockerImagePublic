# coding=utf-8
# 修改CNVkit的Hard-coded参数
# 注意，一般都不用修改
# 除非你有特殊需求

import os
import argparse
import shutil

# 获取当前文件的绝对路径
def get_current_file_path():
    return os.path.dirname(os.path.abspath(__file__))

# 生成 CNVkit 参数文件的 Python 代码
def generate_cnvkit_params_file(file_path, MIN_REF_COVERAGE, MAX_REF_SPREAD, NULL_LOG2_COVERAGE, GC_MIN_FRACTION, GC_MAX_FRACTION, INSERT_SIZE, force_rewrite):
    content = f'''\"\"\"Hard-coded parameters for CNVkit. These should not change between runs.\"\"\"
# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = {MIN_REF_COVERAGE}
MAX_REF_SPREAD = {MAX_REF_SPREAD}
NULL_LOG2_COVERAGE = {NULL_LOG2_COVERAGE}

# Thresholds used in GC-content masking of bad bins at 'fix' step
GC_MIN_FRACTION = {GC_MIN_FRACTION}
GC_MAX_FRACTION = {GC_MAX_FRACTION}

# Fragment size for paired-end reads
INSERT_SIZE = {INSERT_SIZE}

# Target/bin names that are not meaningful gene names
# (In some UCSF panels, "CGH" probes denote selected intergenic regions)
IGNORE_GENE_NAMES = ("-", ".", "CGH")
ANTITARGET_NAME = "Antitarget"
ANTITARGET_ALIASES = (ANTITARGET_NAME, "Background")

# PAR1/2 start/end definitions
PSEUDO_AUTSOMAL_REGIONS = {{
    "grch37": {{"PAR1X": [60000, 2699520], "PAR2X": [154931043, 155260560], "PAR1Y": [10000, 2649520], "PAR2Y": [59034049, 59363566] }},
    "grch38": {{"PAR1X": [10000, 2781479], "PAR2X": [155701382, 156030895], "PAR1Y": [10000, 2781479], "PAR2Y": [56887902, 57217415] }},
}}
SUPPORTED_GENOMES_FOR_PAR_HANDLING = PSEUDO_AUTSOMAL_REGIONS.keys()
'''

    # 已存在文件，是否重写
    if os.path.exists(file_path):
        if not force_rewrite:
            return
    
    with open(file_path, 'w') as file:
        file.write(content)

# reference的自动模型首选从antitarget修改为target
def modify_reference_auto_model(file_path, new_file_path):
    with open(file_path, 'r') as ori, open(new_file_path, 'w') as new:
        for line in ori:
            if "preferring antitargets" in line:
                new.write(line.replace("preferring antitargets", "preferring targets"))
            elif "sexes[sid] = a_is_xx" in line:
                new.write(line.replace("sexes[sid] = a_is_xx", "sexes[sid] = t_is_xx"))
            else:
                new.write(line)
    # 先保存一个备份，不删除原文件
    shutil.copy(file_path, file_path + ".bak")

    # 将旧的替换成新的
    shutil.move(new_file_path, file_path)

# 传参模式
def main():
    parser = argparse.ArgumentParser(description="Generate CNVkit parameters file")
    parser.add_argument("--file_path", type=str, help="Path to the output file", default="/opt/conda/lib/python3.10/site-packages/cnvlib/params.py")
    parser.add_argument("--MIN_REF_COVERAGE", type=str, help="Minimum reference coverage", default="-5.0")
    parser.add_argument("--MAX_REF_SPREAD", type=str, help="Maximum reference spread", default="1.0")
    parser.add_argument("--NULL_LOG2_COVERAGE", type=str, help="NULL log2 coverage", default="-20.0")
    parser.add_argument("--GC_MIN_FRACTION", type=str, help="Minimum GC fraction", default="0.3")
    parser.add_argument("--GC_MAX_FRACTION", type=str, help="Maximum GC fraction", default="0.7")
    parser.add_argument("--INSERT_SIZE", type=str, help="Insert size", default="250")
    parser.add_argument("--force_rewrite", type=bool, help="Force rewrite", default=False)
    parser.add_argument("--reference_auto_model", type=bool, help="Change reference auto model to prefer targets", default=False)

    args = parser.parse_args()
    generate_cnvkit_params_file(args.file_path, args.MIN_REF_COVERAGE, args.MAX_REF_SPREAD, args.NULL_LOG2_COVERAGE, args.GC_MIN_FRACTION, args.GC_MAX_FRACTION, args.INSERT_SIZE, args.force_rewrite)

    if args.reference_auto_model:
        print("Reference auto model is True, modifying reference auto model")
        modify_reference_auto_model("/opt/conda/lib/python3.10/site-packages/cnvlib/reference.py", "/opt/conda/lib/python3.10/site-packages/cnvlib/reference_new.py")

    # 没有输入时参数时打印目前所有的参数
    print("Usage: python cnvkit_params_modify.py --file_path <path_to_output_file> --force_rewrite True\n\n")
    print("Change reference auto model to prefer targets: python cnvkit_params_modify.py --reference_auto_model True")
    print("Current parameters:", f"{args.file_path}")
    print(f"MIN_REF_COVERAGE: {args.MIN_REF_COVERAGE}")
    print(f"MAX_REF_SPREAD: {args.MAX_REF_SPREAD}")
    print(f"NULL_LOG2_COVERAGE: {args.NULL_LOG2_COVERAGE}")
    print(f"GC_MIN_FRACTION: {args.GC_MIN_FRACTION}")
    print(f"GC_MAX_FRACTION: {args.GC_MAX_FRACTION}")
    print(f"INSERT_SIZE: {args.INSERT_SIZE}")

# 运行传参模式
if __name__ == "__main__":
    main()
