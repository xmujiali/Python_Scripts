##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import re
import os
import configparser

parms = configparser.ConfigParser()
parms.read('setting.ini')
total_carbons = parms.getint('topo', 'total_carbons')
# a more eval() to eliminate the extra quotation marks
GAUSSIAN = eval(parms.get('gaussian', 'GAUSSIAN'))

line = ' ---------------------------------------------------------------------\n'
line_1 = r'\s+Input orientation:\s+\n'
line_2 = r'\s+Standard orientation:\s+\n' 

reg_inp = re.compile(line_1 + line + r'.*?' + line + r'(.*?)' + line, re.S)
# reg_std = re.compile(line_2 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_select_coords = re.compile(r'^\s*\d+\s*(\d+)\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)$', re.M)

adjacent_tab = np.zeros((total_carbons+1,3), dtype=np.int32)
adjacent_tab_calc = np.zeros((total_carbons+1,3), dtype=np.int32)
coords = np.zeros((total_carbons+1,3), dtype=float)
dist_mat = np.zeros((total_carbons+1,total_carbons+1), dtype=float)

def check_adj_tab(name):
    # get adj_tab
    with open(f'..{os.sep}coords{os.sep}C_{total_carbons}_{isomer}', 'r') as g:
            
        for line in g.readlines()[total_carbons+2:2*total_carbons+2]:
            center, neighbors = line.split(':')
            center = int(center)
            nbr1,nbr2,nbr3 = neighbors.split(',')
            nbr1 = int(nbr1)
            nbr2 = int(nbr2)
            nbr3 = int(nbr3)
            nbr_list = [nbr1,nbr2,nbr3]
            nbr_list.sort()
            adjacent_tab[center] = np.array(nbr_list)
        

    with open(name + '.log', 'r') as f_log:
        text = f_log.read()
        tmp = re.findall(reg_inp, text)[-1]
        tmp = re.findall(reg_select_coords, tmp)
        for i in range(total_carbons):
            _tmp = tmp[i]
            coords[i+1][0] = float(_tmp[1])
            coords[i+1][1] = float(_tmp[2])
            coords[i+1][2] = float(_tmp[3])
    
    
    for i in range(1,total_carbons):
        for j in range(i,total_carbons+1):
            # looping for i < j, 
            # calculate the distance between i and j
            r_ij = np.linalg.norm(coords[i] - coords[j])
            dist_mat[i][j] = dist_mat[j][i] = r_ij
    
    for i in range(1,total_carbons+1):
        # There are two 'index' must be excluded,
        # namely, 0 and i (itself)
        nbr_list = dist_mat[i].argsort()[2:5]
        nbr_list.sort()
        adjacent_tab_calc[i] = np.array(nbr_list)
    
    # print(name)
    # print('coords = \n', coords)
    # print('dist_mat = \n', dist_mat)
    # print('adj_tab_calc = \n', adjacent_tab_calc)
    # print('adj_tab =\n', adjacent_tab)
    # print(adjacent_tab_calc - adjacent_tab)
    if np.count_nonzero((adjacent_tab_calc - adjacent_tab)) == 0:
        return True
    else:
        return False



ALL_succ = True
with open('isomers', 'r') as f_isomer, \
    open('adj_error', 'w') as f_failed, \
      open('redo_opt.sh', 'w') as f_redo:
    for line in f_isomer.readlines():
        isomer = int(line.strip())
        name = f'C_{total_carbons}_{isomer}'
        if not check_adj_tab(name):
            ALL_succ = False
            f_failed.write(f'{isomer}\n')
            f_redo.write(f'{GAUSSIAN} C_{total_carbons}_{isomer}.gjf')
        
if ALL_succ:
    print('All optimization are normally terminated!')
else:
    print('Warning: isomers incompatible with its adjacent table are listed in file adj_error.')
