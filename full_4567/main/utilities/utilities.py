##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import parms
import numpy as np
import re

def skip_lines(f, n):
    for i in range(n):
        f.readline()

def mkdir_safe(d):
    if os.path.exists(d):
        print(f'{d} is already exists!')
        tmp = input('press d to delete it,\npress c to continue,\nor press any other key to exit!\n')
        tmp.lower()
        if tmp == 'd' or tmp == 'D':
            shutil.rmtree(d)
        elif tmp == 'c' or tmp == 'C':
            return
        else:
            exit(-1)
    os.mkdir(d)

def mkdir_force(d):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.mkdir(d)

def menu(promote_string, lis=None, custom=True):
    print(promote_string)
    if lis:
        for i,d in enumerate(lis):
            print('\t', i+1, '\t', d)
        if custom:
            x = input('or type your custom one.\n')
        else:
            x = input('choose one above:\n')  
        
        if x in [str(i+1)  for i in range(len(lis))]:
            return lis[int(x)-1]
        elif custom:
            return x
        else:
            return menu(promote_string, lis, custom)
    else:
        raise ValueError('No choices in argument lis')

def list_files(dir, extension=''):
    result = []
    for root, dirs, files in os.walk(dir):
        for f in files:
            n = len(extension)
            if (not extension) or (extension and f[-n:]==extension):
                result.append(root + os.sep + f)
    return result

# read adjacent table
def read_adj_tab(isomer, adjacent_tab):
    with open(f'{parms.main_path}{os.sep}{parms.CAGE_data_path}{os.sep}C_{parms.total_carbons}_{isomer}', 'r') as g:
        for line in g.readlines()[parms.total_carbons+2:2*parms.total_carbons+2]:
            center, neighbors = line.split(':')
            center = int(center)
            nbr1,nbr2,nbr3 = neighbors.split(',')
            nbr1 = int(nbr1)
            nbr2 = int(nbr2)
            nbr3 = int(nbr3)
            nbrs = [nbr1,nbr2,nbr3]
            nbrs.sort()
            adjacent_tab[center] = np.array(nbrs)


def write_gjf(isomer, template, coords=None):
    gjf_name = f'C_{parms.total_carbons}_{isomer}.gjf'
    # Needs writing coordinates ?
    coord_flag = re.findall(r'<geometry>', template)

    with open(gjf_name, 'w') as f_gjf:
        result = re.sub('<total_carbons>', str(parms.total_carbons), 
        template, re.S)
        result = re.sub('<isomer>', str(isomer), result, re.S)
        if coord_flag is not None:
            tmp = []
            if coords is None:
                with open(f'{parms.main_path}{os.sep}{parms.CAGE_data_path}{os.sep}C_{parms.total_carbons}_{isomer}', 'r') as f:
                    for line in f.readlines()[:parms.total_carbons]:
                        tmp.append(line)
            else:
                for i in range(1, parms.total_carbons+1):
                    tmp.append(f'C\t{coords[i][0]: 10.6f}\t{coords[i][1]: 10.6f}\t{coords[i][2]: 10.6f}\n')
            geom_string = ''.join(tmp)
            result = re.sub(r'<geometry>', geom_string, result, re.S)
        f_gjf.write(result)




if __name__ == '__main__':
    pass