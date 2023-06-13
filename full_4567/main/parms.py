##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import logging
import configparser

'''
parms introduces and manages all configure options all in one.
users only need import parms and use parms.xxx to access variable xxx
TODO: value check for options
'''



# initialize paramters
parms = configparser.ConfigParser()

if os.path.exists(f'.{os.sep}setting.ini'):
       parms.read(f'.{os.sep}setting.ini', encoding='UTF-8')
elif os.path.exists(f'..{os.sep}setting.ini'):
       parms.read(f'..{os.sep}setting.ini', encoding='UTF-8')
else:
        print('Missing configure file setting.ini!')
        exit(-1) 
       
# topo
MIN_RING_SIZE = parms.getint('topo', 'MIN_RING_SIZE')
MAX_RING_SIZE = parms.getint('topo', 'MAX_RING_SIZE')
total_carbons = parms.getint('topo', 'total_carbons')
total_isomers = parms.getint('topo', 'total_isomers')


# path_and_name
main_path = eval(parms.get('path_and_name', 'main_path'))
CAGE_data_path = eval(parms.get('path_and_name', 'CAGE_data_path'))
w3d = eval(parms.get('path_and_name', 'w3d'))
topology = eval(parms.get('path_and_name', 'topology'))
kekule = eval(parms.get('path_and_name', 'kekule'))
huckel = eval(parms.get('path_and_name', 'huckel'))
directories = eval(parms.get('path_and_name', 'directories'))
gjf_template_path = eval(parms.get('path_and_name', 'gjf_template_path'))
opt_error_log = eval(parms.get('path_and_name', 'opt_error_log'))
adj_error_log = eval(parms.get('path_and_name', 'adj_error_log'))
nlo_error_log = eval(parms.get('path_and_name', 'nlo_error_log'))
result = eval(parms.get('path_and_name', 'result'))
logging_file = eval(parms.get('path_and_name', 'logging_file'))



# w3d
ISOMER_HEAD_LINE = parms.getint('w3d', 'ISOMER_HEAD_LINE')
ISOMER_TAIL_LINE = parms.getint('w3d', 'ISOMER_TAIL_LINE')

# gaussian
GAUSSIAN = eval(parms.get('gaussian', 'GAUSSIAN'))
charge = parms.getint('gaussian', 'charge')
spin_multiplicity = parms.getint('gaussian', 'spin_multiplicity')

# reconstructor
Max_Tries = parms.getint('reconstructor', 'Max_Tries')
Max_Cycle = parms.getint('reconstructor', 'Max_Cycle')
R_nonbond_min = parms.getfloat('reconstructor', 'R_nonbond_min')
R_bond_max = parms.getfloat('reconstructor', 'R_bond_max')
R_bond_min = parms.getfloat('reconstructor', 'R_bond_min')
Factor_change = parms.getfloat('reconstructor', 'Factor_change')


# other
report_frequency = parms.getint('other', 'report_frequency')


# level =
# logging.DEBUG
# logging.INFO
# logging.WARNING
# logging.ERROR
# logging.CRITICAL
logging.basicConfig(filename=main_path+os.sep+logging_file,
                    filemode='a',
                    format='%(asctime)s - %(filename)s - %(message)s',
                    level=logging.DEBUG)
# debug = logging.debug
# info = logging.info
# warning = logging.warning
# error = logging.error
# critical = logging.critical
