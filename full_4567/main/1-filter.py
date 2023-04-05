##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""
import re
import os
import shutil
import parms

print('Choose a filter from following:\n')
print('-'*20)
while True:
    for k,v in parms.filters.items():
        print('key = ', k)
        des = v['description']
        print(des)
        print('-'*20)
    key = input('Which one do you want to use: ')
    if key in parms.filters.keys():
        isomers = parms.filters[key]['filter']()
        break
    else:
        print('Wrong choice. Try again!')

while True:
    template = input('Give me a template file for gaussian inputs:')
    if os.path.exists(template):
        with open(template, 'r') as f:
            contents = f.read()
            break
    else:
        print('File not exits. Try again!')


dir_name = input('Input the directory name where computational jobs located:\n')
if os.path.exists(dir_name):
    shutil.rmtree(dir_name)
os.mkdir(dir_name)
os.chdir(dir_name)

with open('isomers', 'w') as f_iso, open('run_jobs.sh', 'w') as f_run:
    # if there is NO %oldchk line, then we need write coordinates
    coord_flag = re.match('%oldchk', contents, re.S)
    if coord_flag is None:
        coord_path = input('Give me the location of coodinate files[enter for ../coords]:')
        if coord_path == '':
            coord_path = '..' + os.sep + 'coords'

    for i in isomers:
        f_iso.write(str(i) + '\n')
        name = 'C_%d_%s' % (parms.total_carbons, i)
        gjf_name = name + '.gjf'
        f_run.write(parms.GAUSSIAN + ' ' + gjf_name + '\n')

        with open(gjf_name, 'w') as f_gjf:
            result = re.sub("\(isomer\)", str(i), contents, re.S)
            f_gjf.write(result)

            # if there is NO %oldchk line, then we need write coordinates
            if coord_flag is None:
                coord_name = coord_path + os.sep + name
                with open(coord_name, 'r') as coord:
                    for line in coord.readlines()[:parms.total_carbons]:
                        f_gjf.write(line)
            f_gjf.write('\n')

print('Job done!')