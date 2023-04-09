##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""
import pandas as pd
import os
import re
import parms

result_file = input('Give me the file name to store optimized energies[Enter for ../result.csv]:')
if result_file == '':
    result_file = '..' + os.sep + 'result.csv'
if os.path.exists(result_file):
    res = pd.read_csv(result_file, index_col=0)
else:
    res = pd.DataFrame()

col_name = input('column name = ')
res[col_name] = 999999999.0


with open('isomers', 'r') as f_iso:
    for line in f_iso.readlines():
        isomer = line[:-1]
        name = 'C_%d_%s' % (parms.total_carbons, isomer)
        with open(name + '.log', 'r') as f_log:
            energies = []
            for line in f_log.readlines():
                if re.match('^ SCF Done:', line):
                    energy = float(line.split()[4])
                    energies.append(energy)
            res.loc[int(isomer), col_name] = energies[-1]
res.to_csv(result_file)
