##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import parms
from utilities.reconstructor import reconstruct 

def check_opt(isomer):
    log_name = f'C_{parms.total_carbons}_{isomer}.log'
    if os.path.exists(log_name):
        with open(log_name, 'r') as log:
            content = log.read()
        if re.findall(r' -- Stationary point found.', content):
            return True
    return False

if __name__ == '__main__':
    ALL_succ = True
    with open('isomers', 'r') as fi, open(parms.opt_error_log, 'w') as feol:
        for line in fi.readlines():
            isomer = int(line.strip())
            if isomer % parms.report_frequency == 1:
                print(f'Handling isomer {isomer}')
            if not check_opt(isomer):
                feol.write(f'{isomer}\n')
                ALL_succ = False
    if ALL_succ:
        os.remove(parms.opt_error_log)
        print('All optimization are normally terminated!')
    else:
        print(f'FAILED optimization jobs are logged in {parms.opt_error_log}.')

