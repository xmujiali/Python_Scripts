##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import re
import os
import parms
import utilities.utilities as util

line = ' ---------------------------------------------------------------------\n'
line_1 = r'\s+Input orientation:\s+\n'
line_2 = r'\s+Standard orientation:\s+\n' 

reg_inp = re.compile(line_1 + line + r'.*?' + line + r'(.*?)' + line, re.S)
# reg_std = re.compile(line_2 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_select_coords = re.compile(r'^\s*\d+\s*(\d+)\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)$', re.M)

adjacent_tab = np.zeros((parms.total_carbons+1,3), dtype=np.int32)
adjacent_tab_calc = np.zeros((parms.total_carbons+1,3), dtype=np.int32)
coords = np.zeros((parms.total_carbons+1,3), dtype=float)
dist_mat = np.zeros((parms.total_carbons+1,parms.total_carbons+1), dtype=float)

def adj_check(isomer):
    util.read_adj_tab(isomer, adjacent_tab)
        
    log_name = f'C_{parms.total_carbons}_{isomer}.log'
    with open(log_name, 'r') as f:
        text = f.read()
        tmp = re.findall(reg_inp, text)[-1]
        tmp = re.findall(reg_select_coords, tmp)
        for i in range(parms.total_carbons):
            _tmp = tmp[i]
            coords[i+1][0] = float(_tmp[1])
            coords[i+1][1] = float(_tmp[2])
            coords[i+1][2] = float(_tmp[3])

    for i in range(1,parms.total_carbons):
        for j in range(i,parms.total_carbons+1):
            # looping for i < j, 
            # calculate the distance between i and j
            r_ij = np.linalg.norm(coords[i] - coords[j])
            dist_mat[i][j] = dist_mat[j][i] = r_ij
    
    for i in range(1,parms.total_carbons+1):
        # There are two 'index' must be excluded,
        # namely, 0 and i (itself)
        nbr_list = dist_mat[i].argsort()[2:5]
        nbr_list.sort()
        adjacent_tab_calc[i] = np.array(nbr_list)
    
    if 0 == np.count_nonzero((adjacent_tab_calc - adjacent_tab)):
        return True
    else:
        return False


if __name__ == '__main__':
    ALL_succ = True
    with open('isomers', 'r') as fi, open(parms.adj_error_log, 'w') as fael:
        for line in fi.readlines():
            isomer = int(line.strip())
            if isomer % parms.report_frequency == 1:
                print(f'Handling isomer {isomer}')
            if not adj_check(isomer):
                ALL_succ = False
                fael.write(f'{isomer}\n')            
    if ALL_succ:
        os.remove(parms.adj_error_log)
        print('Optimized geometry and adjacent table are compatible for all isomers!')
    else:
        print(f'Warning: isomers with incompatible optimized geometry and adjacent table are listed in file {parms.adj_error_log}.')
