# -*- coding: utf-8 -*-
"""
version: python 3.0
usage: python get_ASprofile_ref_hdrs.py path/species.genome.fa species
"""

import sys
import re
import fileinput
import pandas
import os.path

if len(sys.argv) < 3:
    sys.exit("python error")

FA = sys.argv[1]
# species = os.path.basename(FA).split(".")[0]
species = sys.argv[2]

dic_chr = {}
temp_chr = ''

for line in fileinput.input(FA):
    line = line.strip()
    
    # re.compile函数根据包含的正则表达式的字符串创建模式对象。可以实现更有效率的匹配。
    pat = re.compile(r'^>')
    
    # match()从首字母开始开始匹配，string如果包含pattern子串，则匹配成功，返回Match对象；
    # 失败则返回None，若要完全匹配，pattern要以$结尾。
    match = pat.match(line)
    if match:
        a = line.split(" ")[0]
        temp_chr = a
        
        dic_chr.setdefault(temp_chr, []).append(0)
        dic_chr.setdefault(temp_chr, []).append(0)
    else:
        dic_chr[temp_chr][0] += len(line)
        dic_chr[temp_chr][1] += (line.count("N") + line.count("n"))
    
fileinput.close()


import pandas as pd
import numpy as np

df = pd.DataFrame.from_dict(dic_chr, orient='index')
df = df.reset_index()

for i in range(0, len(df)):
    df.iloc[i, 1] = "/len=" + str(df.iloc[i, 1])
    df.iloc[i, 2] = "/nonNlen=" + str(df.iloc[i, 2])
    
df["species"] = "/org=" + species

np.savetxt(species + '.fa.hdrs', df.values, fmt='%s')