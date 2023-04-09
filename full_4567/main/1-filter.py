##!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import os
import shutil
import configparser
from util_filters import filters

# initialize paramters
parms = configparser.ConfigParser()
parms.read('.\setting.ini', encoding='UTF-8')

MIN_RING_SIZE = parms.getint('topo', 'MIN_RING_SIZE')
MAX_RING_SIZE = parms.getint('topo', 'MAX_RING_SIZE')
total_carbons = parms.getint('topo', 'total_carbons')
total_isomers = parms.getint('topo', 'total_isomers')
ISOMER_HEAD_LINE = parms.getint('w3d', 'ISOMER_HEAD_LINE')
ISOMER_TAIL_LINE = parms.getint('w3d', 'ISOMER_TAIL_LINE')

working_directory = eval(parms.get('name', 'working_directory'))
if os.path.exists(working_directory):
    shutil.rmtree(working_directory)
os.mkdir(working_directory)
os.chdir(working_directory)
sample_gjf_file = eval(parms.get('name', 'sample_gjf_file'))
with open(sample_gjf_file, 'r') as f:
            contents = f.read()
coord_path = eval(parms.get('name', 'coord_path'))


# a more eval() to eliminate the extra quotation marks
GAUSSIAN = eval(parms.get('gaussian', 'GAUSSIAN'))


print('Choose a filter from following:\n')
print('-'*20)
while True:
    for k,v in filters.items():
        print('key = ', k)
        des = v['description']
        print(des)
        print('-'*20)
    key = input('Which one do you want to use: ')
    if key in filters.keys():
        isomers = filters[key]['filter']()
        break
    else:
        print('Wrong choice. Try again!')

with open('isomers', 'w') as f_iso, open('run_jobs.sh', 'w') as f_run:
    # if there is NO %oldchk line, then we need write coordinates
    coord_flag = re.match('%oldchk', contents, re.S)

    for i in isomers:
        f_iso.write(str(i) + '\n')
        name = 'C_%d_%s' % (total_carbons, i)
        gjf_name = name + '.gjf'
        f_run.write(GAUSSIAN + ' ' + gjf_name + '\n')

        with open(gjf_name, 'w') as f_gjf:
            result = re.sub("\(isomer\)", str(i), contents, re.S)
            f_gjf.write(result)

            # if there is NO %oldchk line, then we need write coordinates
            if coord_flag is None:
                coord_name = coord_path + os.sep + name
                with open(coord_name, 'r') as coord:
                    for line in coord.readlines()[:total_carbons]:
                        f_gjf.write(line)
            f_gjf.write('\n')

print('Job done!')