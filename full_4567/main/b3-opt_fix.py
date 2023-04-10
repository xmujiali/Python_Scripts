##!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import parms
import utilities.utilities as util
from utilities.reconstructor import reconstruct 

if __name__ == '__main__':
    ps = '''How to do opt_fix:
    1. by Modifying Routine Section (MRS)
    2. by Reconstruction Geometry from Adjacent Table(RGAT)'''

    response = util.menu('How to do opt_fix:', ['MRS','RGAT'], custom=False)

    if 'MRS' == response:
        reg_exp = '^\s*#\s*[Pp]'
        routine_section = input('Input the new routine section:\n') + '\n'
        with open(parms.opt_error_log,'r') as foel, \
            open('run_opt_fix.sh', 'w') as frof:
            for line in foel.readlines():
                isomer = int(line.strip())
                frof.write(f'{parms.GAUSSIAN} C_{parms.total_carbons}_{isomer}.gjf\n')
                print(f'Handling Isomer: {isomer}')
                gjf_name = f'C_{parms.total_carbons}_{isomer}.gjf'
                with open(gjf_name, 'r') as gjf:
                    lines = gjf.readlines()
                for i,x in enumerate(lines):
                    if re.match(reg_exp, x):
                        lines[i] = routine_section
                with open(gjf_name, 'w') as gjf:
                    gjf.writelines(lines)
                    gjf.write('\n')
    elif 'RGAT' == response:
        with open(parms.opt_error_log,'r') as foel, \
            open('run_opt_fix.sh', 'w') as frof:
            for line in foel.readlines():
                isomer = int(line.strip())
                frof.write(f'{parms.GAUSSIAN} C_{parms.total_carbons}_{isomer}.gjf\n')
                print(f'Handling Isomer: {isomer}')
                reconstruct(isomer)
    print('bash "run_opt_fix.sh" to do fix optimization.')

