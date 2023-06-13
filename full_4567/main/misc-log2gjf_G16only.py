##!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Sep  4 19:04:00 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import numpy as np
import shutil
import os
import re

prefix_in = 'C_28_'
suffix_in = '.log'
prefix_out = 'C_28_'
suffix_out = '.gjf'

template = input('Give me the template file for gaussian inputs:\n')
with open(template, 'r') as te:
    contents = te.read()

job_file = input('Give me the file containing isomers u want to transform form log 2gjf:\n')

gjf_path = 'log2gjf'
if os.path.exists(gjf_path):
    shutil.rmtree(gjf_path)
os.mkdir(gjf_path)


line = ' ---------------------------------------------------------------------\n'
line_1 = r'\s+Input orientation:\s+\n'
line_2 = r'\s+Standard orientation:\s+\n' 

reg_inp = re.compile(line_1 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_std = re.compile(line_2 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_select_coords = re.compile(r'^\s*\d+\s*(\d+)\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)$', re.M)

with open(job_file, 'r') as f:
    for line in f.readlines():
        isomer = line[:-1]
        name_in = prefix_in +  str(isomer) + suffix_in
        name_out = gjf_path + os.sep +  prefix_out +  str(isomer) + suffix_out
        with open(name_in, 'r') as f_in, open(name_out, 'w') as f_out:
            text = f_in.read()
            tmp = re.findall(reg_inp, text)[-1]
            nos = np.zeros(parms.total_carbons, dtype=np.int32)


            tmp = re.findall(reg_select_coords, tmp)
            result = re.sub("\(isomer\)", isomer, contents, re.S)
            f_out.write(result)
            f_out.write('\n')
            for x in tmp:
                f_out.write('%3s,%10s,%10s,%10s\n' % (x[0], x[1], x[2], x[3]))
            f_out.write('\n')
