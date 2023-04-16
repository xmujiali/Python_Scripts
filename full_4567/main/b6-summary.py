##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import re
import logging
import parms

if os.path.exists(parms.result):
    res = pd.read_csv(parms.result, index_col=0)
else:
    res = pd.DataFrame()

col_name = input('column name = ')
res[col_name] = np.nan
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
            res.loc[isomer, col_name] = energies[-1]
res.to_csv(parms.result)
