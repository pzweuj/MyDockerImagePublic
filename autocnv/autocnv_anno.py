# coding=utf-8
# pzw
# 20241224

import sys
import argparse

# 环境变量
sys.path.insert(0, '/opt/autopvs1/autopvs1')
sys.path.insert(0, '/opt/autocnv/autocnv')
from autocnv.annotate import AnnotateHelper

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='CNV注释工具', formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='''
        使用示例:
            python autocnv_anno.py -i input.bed -o output.xlsx

        输入文件格式要求:
            BED格式文件，包含以下列，必须包含标题行，标题行不要以#开头：
            - chr    染色体
            - start  起始位置
            - end    终止位置
            - type   CNV类型 (dup/del)
            - error  错误率(可以不包含)
    '''
    )
    parser.add_argument('-i', '--input', required=True, help='输入bed文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出excel文件路径')
    
    # 如果没有参数，打印帮助信息并退出
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    # 解析参数
    args = parser.parse_args()
    
    # 执行注释
    annotate = AnnotateHelper()
    annotation = annotate.annotation_file(args.input, args.output)

if __name__ == '__main__':
    main()

