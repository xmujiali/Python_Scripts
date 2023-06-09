##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import parms
from utilities.reconstructor import reconstruct


if __name__ == '__main__':
    with open(parms.adj_error_log,'r') as fael, open('run_adj_fix.sh', 'w') as fraf:
        for line in fael.readlines():
            isomer = int(line.strip())
            print(isomer)
            fraf.write(f'{parms.GAUSSIAN} C_{parms.total_carbons}_{isomer}.gjf\n')
            print(f'Handling Isomer: {isomer}')
            reconstruct(isomer)
    logging.info('bash "run_adj_fix.sh" to do fix optimization.')
