##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 01:26:16 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import numpy as np
import pandas as pd
import os
import re
import parms

import logging
#设置日志输出格式
logging.basicConfig(
    #设置日志输出等级
    level=logging.DEBUG,
    # level=logging.WARNING,
    filename="DEBUG.log", # 日志输出的文件位置和文件名
    # 文件的写入格式，w为重新写入文件，默认是追加
    filemode="w",
    # 日志输出的格式
    format="%(asctime)s - %(name)s - %(levelname)-9s - %(filename)-8s : %(lineno)-4s line - %(message)s", 
    # 时间输出的格式
    datefmt="%Y-%m-%d %H:%M:%S"
    )

reg_total = r'\\Dipole=(.*?)\\Polar=(.*?)\\HyperPolar=(.*?)\\'
reg_total = re.compile(reg_total)    


result_file = input('Give me the file name to store NLO coefficients[Enter for ../result.csv]:')
prefix = input('Give me a prefix for computational method, to label NLO data:')
if result_file == '':
    result_file = '..' + os.sep + 'result.csv'
if os.path.exists(result_file):
    res = pd.read_csv(result_file, index_col=0)
else:
    res = pd.DataFrame()


ALL_succ = True
with open('isomers', 'r') as f_iso, open('failed_jobs', 'w') as f_fail, open('redo_jobs.sh', 'w') as f_re:
    for line in f_iso.readlines():
        isomer = line[:-1]
        name = 'C_%d_%s' %(parms.total_carbons, isomer)
        logging.debug('-'*10 + 'handling isomer #'+isomer + '-'*10)
        with open(name + '.log') as log:
            lines = log.readlines()[-500:]
            # Normal termination
            last_line = lines[-1]
            if not re.match('^ Normal termination of Gaussian', last_line):
                ALL_succ = False
                f_fail.write(name + '\n')
                f_re.write((parms.GAUSSIAN + ' ' + name + '.gjf\n'))
            
            # get NLO properties
            text = ''
            for line in lines:
                text += line.strip()
            text = reg_total.search(text)
            
            # handle dipole
            mu_x, mu_y, mu_z = eval(str(text.group(1)))
            mu_vec = np.array([mu_x, mu_y, mu_z])
            mu_total = np.linalg.norm(mu_vec)
            logging.debug('mu_x, mu_y, mu_z, mu_total\n = %f, %f, %f, %f' % (mu_x, mu_y, mu_z, mu_total))

            # handle polar
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

            
            # handle beta
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

            col_name = prefix+'beta'
            res.loc[int(isomer), col_name] = beta_magnitude

            logging.debug('-'*10 + 'Done isomer #' + isomer + '-'*10)
res.to_csv(result_file)
if ALL_succ:
    print('All NLO calculations are normally terminated!\n')
else:
    print('FAILED semi-NLO jobs are logged in failed_jobs.')
    print('Try to fixed them, and run "bash redo_jobs.sh" to restart failed optimizations')

