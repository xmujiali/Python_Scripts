##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
import configparser


def check_opt(name):
    with open(name + '.log') as log:
        log_text = log.readlines()
    for line in log_text:
        if line.strip() == '-- Stationary point found.':
            return True
    return False


parms = configparser.ConfigParser()
parms.read('setting.ini')
total_carbons = parms.getint('topo', 'total_carbons')
# a more eval() to eliminate the extra quotation marks
GAUSSIAN = eval(parms.get('gaussian', 'GAUSSIAN'))


ALL_succ = True
with open('isomers', 'r') as f_isomer, \
    open('opt_err', 'w') as f_failed, \
    open('redo_opt.sh', 'w') as f_redo:

    for line in f_isomer.readlines():
        isomer = int(line.strip())
        name = f'C_{total_carbons}_{isomer}'
        if not check_opt(name):
            f_failed.write(f'{isomer}\n' + '\n')
            f_redo.write(f'{GAUSSIAN} {name}.gjf\n')
            ALL_succ = False
        
        
if ALL_succ:
    print('All optimization are normally terminated!')
else:
    print('FAILED optimization jobs are logged in failed_jobs.')
    print('Try to fixed them, and run "bash redo_opt.sh" to restart failed optimizations')
