# coding=utf-8
# pzw
# 20241223

import os
import sys

# 切换工作目录
sys.path.insert(0, '/opt/autopvs1/autopvs1')
sys.path.insert(0, '/opt/autocnv/autocnv')
from autoPVS1 import AutoPVS1

# 测试
demo = AutoPVS1('13-113803407-G-A', 'GRCh37')
if demo.islof:
    print(demo.hgvs_c, demo.hgvs_p, demo.consequence, demo.pvs1.criterion, 
          demo.pvs1.strength_raw, demo.pvs1.strength)
    
