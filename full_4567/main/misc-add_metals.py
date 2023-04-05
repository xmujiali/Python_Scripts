##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import os
import sys
import numpy as np
import parms


def skip_lines(f: object, n: int):
    for i in range(n):
        f.readline()


def make_run1(name, coords, **kw):
    with open(name + '.gjf', 'w') as gjf:
        gjf.write('%chk=' + name + '.chk\n')
        gjf.write('%%mem=%s\n' % kw['mem'])
        gjf.write('%%nproc=%s\n' % kw['nproc'])
        gjf.write('#T %s opt=%s scf=%s nosymm\n'
                  % (kw['method'], kw['opt'], kw['scf']))
        gjf.write('\n')
        gjf.write(name + '\n')
        gjf.write('\n')
        spin_multiplicity = parms.charge % 2 + 1
        gjf.write('%d, %d\n' % (parms.charge, spin_multiplicity))
        gjf.writelines(coords)
        gjf.write('\n')

##############################################################################

# allocate some arrays to store information for each isomer
coords = np.zeros((parms.total_carbons+1, 3))
adj_tab = np.zeros((parms.total_carbons+1, 3), dtype='int32')

print('coords location:[../coords]')
dir_coords = input()
print('DFT results location:[../DFT_calc]')
dir_dft = input()
isomer = int(input('Input the #isomer you want to handle:'))

os.path.join(dir_coords, '')

with open(os.path, 'r') as f_c, open(
        sys.argv[2], 'r') as f_l, open('out.gjf', 'w') as f_o:
    # write the head part of gjf file
    f_o.write('%chk=xxx.chk\n')
    f_o.write('%%mem=10000MB\n')
    f_o.write('%nproc=8\n')
    f_o.write('#T M062X/6-31G* opt=verytight scf=(xqc,tight) nosymm\n')
    f_o.write('\n')
    f_o.write('AAA\n')
    f_o.write('\n')
    f_o.write('0\t1\n')


    # read coordinates fromm log file
    flag1 = False
    flag2 = False
    pass_lines = 5
    read_lines = parms.total_carbons
    for line in f_l.readlines():
        if line.strip() == '-- Stationary point found.':
            flag1 = True
        if flag1 and line.strip() == 'Standard orientation:':
            flag2 = True
        if flag1 and flag2:
            if pass_lines > 0:
                pass_lines -= 1
            elif read_lines > 0:
                read_lines -= 1
                temp = line.split()
                if 6 == int(temp[1]):
                    car = int(temp[0])
                    coords[car][0] = float(temp[3])
                    coords[car][1] = float(temp[4])
                    coords[car][2] = float(temp[5])
                    f_o.write('C\t%9.6f\t%9.6f\t%9.6f\n'
                              % tuple(coords[car]))

    # read the adj_tab
    skip_lines(f_c, parms.total_carbons + 4)
    for i in range(parms.total_carbons):
        temp = f_c.readline().split(':')
        car = int(temp[0])
        nbrs = temp[1].split(',')
        adj_tab[car][0] = int(nbrs[0])
        adj_tab[car][1] = int(nbrs[1])
        adj_tab[car][2] = int(nbrs[2])

    # read the bond2type
    skip_lines(f_c, 14)
    bond2type = eval(f_c.readline())
    for c1,c2 in bond2type.keys():
        print(c1,c2)




    # write the tail part of gjf file
    f_o.write('\n')
    f_o.write('Y\t0\n')
    f_o.write('SDD\n')
    f_o.write('****\n')
    f_o.write('C\t0\n')
    f_o.write('6-31G*\n')
    f_o.write('****\n')
    f_o.write('\n')
    f_o.write('Y\t0\n')
    f_o.write('SDD\n')

