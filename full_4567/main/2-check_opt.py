##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import parms


def check_opt(name):
    success = False
    with open(name + '.log') as log:
        for line in log.readlines():
            if line.strip() == '-- Stationary point found.':
                success = True
                break
    if success:
        return True
    else:
        return False


ALL_succ = True
with open('isomers', 'r') as f_isomer, open('failed_jobs', 'w') as f_failed, open('redo_jobs.sh', 'w') as f_redo:
    for line in f_isomer.readlines():
        isomer = line[:-1]
        name = 'C_%d_%s' % (parms.total_carbons, isomer)
        if not check_opt(name):
            f_failed.write(line)
            name = 'C_%d_%s' % (parms.total_carbons, isomer)
            f_redo.write(parms.GAUSSIAN + ' ' + name + '.gjf\n')
            ALL_succ = False

if ALL_succ:
    print('All optimization are normally terminated!')
else:
    print('FAILED optimization jobs are logged in failed_jobs.')
    print('Try to fixed them, and run "bash redo_jobs.sh" to restart failed optimizations')
