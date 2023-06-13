##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import re
import logging
import parms

if __name__ == '__main__':
    file = parms.main_path + os.sep + parms.result
    if os.path.exists(file):
        res = pd.read_csv(file, index_col=0)
    else:
        res = pd.DataFrame()

    pre_fix = input('pre_fix of column name:')
    col_name1 = pre_fix + '_au'
    res[col_name1] = int(input('default value = ')) 
    with open('isomers', 'r') as f_iso:
        for line in f_iso.readlines():
            isomer = int(line.strip())

            if 1 == isomer % parms.report_frequency:
                logging.info(f'Handling isomer {isomer}')

            log_name = f'C_{parms.total_carbons}_{isomer}.log'
            with open(log_name, 'r') as f:
                energies = []
                for line in f.readlines():
                    if re.match('^ SCF Done:', line):
                        energy = float(line.split()[4])
                        energies.append(energy)
                res.loc[isomer, col_name1] = energies[-1]

    E_min = res[col_name1].min()
    col_name2 = pre_fix + '_kcal'
    res[col_name2] = (res[col_name1] - E_min) * 627.51

    res.to_csv(file)
