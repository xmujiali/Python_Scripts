##!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import os
import re
import logging
import parms


if os.path.exists(parms.result):
    res = pd.read_csv(parms.result, index_col=0)
else:
    res = pd.DataFrame()

pre_fix = input('column name =\n')
col_beta = pre_fix+'beta'

# Normal termination of Gaussian xx
reg_NToG = re.compile('^ Normal termination of Gaussian')
# Dipole_Polar_HyperPolar
reg = re.compile(r'\\Dipole=(.*?)\\Polar=(.*?)\\HyperPolar=(.*?)\\')

def nlo_to_beta(lines):
    '''
    analyze the last few lines (summary part) of gaussian .log file,
    to obtain NLO properties mu, alpha and beta. 
    '''
    # if Normal termination
    last_line = lines[-1]
    if not re.match(reg_NToG, last_line):
        return False
    
    # get NLO properties
    text = ''.join(map(str.strip, lines))
    text = reg.search(text)
    if text is None:
        return False
    
    # handling dipole
    mu_x, mu_y, mu_z = eval(str(text.group(1)))
    mu_vec = np.array([mu_x, mu_y, mu_z])
    mu_total = np.linalg.norm(mu_vec)
    logging.debug('mu_x, mu_y, mu_z, mu_total\n = %f, %f, %f, %f' % (mu_x, mu_y, mu_z, mu_total))

    # handling polar
    a_xx, a_xy, a_yy, a_xz, a_yz, a_zz = eval(str(text.group(2)))
    a = np.array([[a_xx, a_xy, a_xz], [a_xy, a_yy, a_yz], [a_xz, a_yz, a_zz]])
    a_iso = np.trace(a)/3
    temp1 = (a_xx-a_yy)**2+(a_xx-a_zz)**2+(a_yy-a_zz)**2
    temp2 = 6*(a_xy**2+a_yz**2+a_xz**2)
    a_aniso1 = np.sqrt((temp1+temp2)/2)
    a_aniso2 = np.sqrt(temp1/2)
    ev = np.linalg.eigvals(a)
    ev.sort()
    a_aniso3 = ev[2] - (ev[0]+ev[1])/2
    logging.debug('a_xx, a_xy, a_yy, a_xz, a_yz, a_zz\n = %f, %f, %f, %f, %f, %f' % (a_xx, a_xy, a_yy, a_xz, a_yz, a_zz))
    logging.debug('a_iso, a_aniso1, a_aniso2, a_aniso3\n = %f, %f, %f, %f' % (a_iso, a_aniso1, a_aniso2, a_aniso3))        

    
    # handling beta
    b_xxx, b_xxy, b_xyy, b_yyy, b_xxz, b_xyz, b_yyz, b_xzz, b_yzz, b_zzz = eval(str(text.group(3)))
    beta_x = -(b_xxx + b_xyy + b_xzz)
    beta_y = -(b_xxy + b_yyy + b_yzz)
    beta_z = -(b_xxz + b_yyz + b_zzz)
    beta_vec = np.array([beta_x, beta_y, beta_z])
    beta_magnitude = np.linalg.norm(beta_vec)
    
    if mu_total > 10**(-6):
        beta_para_mu = np.dot(beta_vec, mu_vec)/mu_total
    else:
        beta_para_mu = np.nan
    
    beta_para_z = 0.6 * beta_z
    beta_orth_z = 0.2 * beta_z
    # beta_orthogonal_z = 0.2*(b_xxz + b_yyz + b_zzz)

    
    logging.debug('b_xxx, b_xxy, b_xyy, b_yyy, b_xxz, b_xyz, b_yyz, b_xzz, b_yzz, b_zzz\n = %f, %f, %f, %f, %f, %f, %f, %f, %f, %f' % (b_xxx, b_xxy, b_xyy, b_yyy, b_xxz, b_xyz, b_yyz, b_xzz, b_yzz, b_zzz))
    logging.debug('beta_x, beta_y, beta_z, beta_magnitude\n = %f, %f, %f, %f' % (beta_x, beta_y, beta_z, beta_magnitude))
    return beta_magnitude





if __name__ == '__main__':
    ALL_succ = True
    with open('isomers', 'r') as fi, \
        open(parms.nlo_error_log, 'w') as fnel, \
            open('run_nlo_fix.sh', 'w') as frnf:
        for line in fi.readlines():
            isomer = int(line.strip())

            log_name = f'C_{parms.total_carbons}_{isomer}.log'
            
            if 1 == isomer % parms.report_frequency:
                logging.info(f'Handling isomer {isomer}')

            succ = True
            with open(log_name, 'r') as log:
                lines = log.readlines()[-500:]

            tmp = nlo_to_beta(lines)
            if tmp:
                ALL_succ = False
                fnel.write(f'{isomer}\n')
                frnf.write(f'{parms.GAUSSIAN} {isomer}.gjf\n')
            else:
                res.loc[isomer, col_beta] = tmp
    
    res.to_csv(parms.result)
    if ALL_succ:
        logging.info('All NLO calculations are normally terminated!')
    else:
        logging.info(f'FAILED semi-NLO jobs are logged in {parms.nlo_error_log}.')
        logging.info('Try to fixed them, and run "bash run_nlo_fix.sh" to restart failed optimizations')

