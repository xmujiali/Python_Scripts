##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""
import numpy as np
import pandas as pd
import os
import re
import parms


line = ' ---------------------------------------------------------------------\n'
line_1 = r'\s+Input orientation:\s+\n'
line_2 = r'\s+Standard orientation:\s+\n' 

reg_inp = re.compile(line_1 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_std = re.compile(line_2 + line + r'.*?' + line + r'(.*?)' + line, re.S)
reg_select_coords = re.compile(r'^\s*(\d+)\s*(\d+)\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)$', re.M)




# 在I1~I4之间的文本是input orientation下，mu, alpha和beta对应的输出文本
# 在D1~D4之间的文本是dipole orientation下，mu, alpha和beta对应的输出文本
I1 = ' Electric dipole moment \(input orientation\):\n'
I2 = ' Dipole polarizability, Alpha \(input orientation\).\n'
I3 = ' First dipole hyperpolarizability, Beta \(input orientation\).\n'
I4 = ' ----------------------------------------------------------------------\n'
D1 = ' Electric dipole moment \(dipole orientation\):\n'
D2 = ' Dipole polarizability, Alpha \(dipole orientation\).\n'
D3 = ' First dipole hyperpolarizability, Beta \(dipole orientation\).\n'
D4 = ' ----------------------------------------------------------------------\n'
# 获取各段文本的正则表达式
reg_partition = re.compile(I1 + r'(.*?)' + I2 + r'(.*?)' + I3 + r'(.*?)' + I4 + r'.*?' + D1 + r'(.*?)' + D2 + r'(.*?)' + D3 + r'(.*?)' + D4, re.S)
# 从每段文本中得到含有数据的行，并且提取该数据的名称与数值
reg_select_nlos = re.compile(r'^\s*(\S*)\s*(\S+D\S*)\s*\S+D\S*\s*\S+D\S*$', re.M)

def get_coeff(d, text):
    tmp = re.findall(reg_select_nlos, text)
    for k,v in tmp:
        d[k] = float(v.replace('D', 'E'))

def get_coords(nos, atnos, coords, text):
    tmp = re.findall(reg_select_coords, text)
    for x in tmp:
        no = int(x[0])
        nos[no-1] = no
        atnos[no-1] = int(x[1])
        coords[no-1] = (float(x[2]), float(x[3]), float(x[4]), 1.0)




###################################################################################

result_file = input('Give me the file name to store NLO coefficients[Enter for ../result.csv]:')
if result_file == '':
    result_file = '..' + os.sep + 'result.csv'
if os.path.exists(result_file):
    res = pd.read_csv(result_file, index_col=0)
else:
    res = pd.DataFrame()

suffix = input('Give me a suffix for computational method, to label NLO data:')
col_name = 'beta_mag_' + suffix
res[col_name] = 0.0


with open('isomers', 'r') as f_iso:
    for line in f_iso.readlines():
        isomer = line[:-1]
        name = 'C_%d_%s' % (parms.total_carbons, isomer)

        with open(name + '.log', 'r') as f_log:
            text = f_log.read()
            tmp = re.findall(reg_inp, text)[0]
            nos = np.zeros(parms.total_carbons, dtype=np.int32)
            atnos = np.zeros(parms.total_carbons, dtype=np.int32)
            coords_inp = np.zeros((parms.total_carbons, 4))
            get_coords(nos, atnos, coords_inp, tmp)


            tmp = re.findall(reg_std, text)[0]
            nos = np.zeros(parms.total_carbons, dtype=np.int32)
            atnos = np.zeros(parms.total_carbons, dtype=np.int32)
            coords_std = np.zeros((parms.total_carbons, 4))
            get_coords(nos, atnos, coords_std, tmp)
            
            
            pinv_inp = np.linalg.pinv(coords_inp)
            A =  pinv_inp @ coords_std 
            result = coords_inp @ A

            # 旋转矩阵
            R = A[0:3, 0:3]
            # print('R = \n', R)
            # print('RT @ R = \n', R.transpose() @ R)
            # print('R @ RT = \n', R @ R.transpose())
            # print('R^-1 @ R \n', np.linalg.inv(R) @ R )
            # print('R @ R^-1 \n', R @ np.linalg.inv(R))

            tmp = re.findall(reg_partition, text)
            if tmp:
                res.loc[int(isomer), col_name] = np.linalg.norm(0)
                print('Isomer #%s has no nlo data found, use zero instead!')
            else:
                # 下面是六段需要的文本,第0段是input orientation下的mu，其余类推
                mu_inp_text, alpha_inp_text, beta_inp_text, mu_dipole_text, alpha_dipole_text, beta_dipole_text = tmp[0]
                # 下面六个字典用来存放读取的数值,第0段是input orientation下的mu，其余类推
                mu_inp, alpha_inp, beta_inp, mu_dipole, alpha_dipole, beta_dipole = {},{},{},{},{},{}

                get_coeff(mu_inp, mu_inp_text)
                get_coeff(alpha_inp, alpha_inp_text)
                get_coeff(beta_inp, beta_inp_text)
                get_coeff(mu_dipole, mu_dipole_text)
                get_coeff(alpha_dipole, alpha_dipole_text)
                get_coeff(beta_dipole, beta_dipole_text)


                # 以下检验我们所得数据与MultiWFN是否一致
                # mu
                mu1 = np.array([mu_inp['x'], mu_inp['y'], mu_inp['z']])
                mu2 = mu1 @ R
                print('mu2 = \n', mu2)

                # alpha
                alpha1 = np.zeros((3,3))
                alpha1[0,0] = alpha_inp['xx']
                alpha1[1,1] = alpha_inp['yy']
                alpha1[2,2] = alpha_inp['zz']
                alpha1[0,1] = alpha1[1,0] = alpha_inp['yx']
                alpha1[0,2] = alpha1[2,0] = alpha_inp['zx']
                alpha1[1,2] = alpha1[2,1] = alpha_inp['zy']
                alpha2 = np.einsum('ij,ik,jl->kl', alpha1, R, R)
                print('alpha3 = \n', alpha2)


                # beta
                beta1 = np.zeros((3,3,3))
                beta1[0,0,0] = beta_inp['xxx']
                beta1[1,1,1] = beta_inp['yyy']
                beta1[2,2,2] = beta_inp['zzz']
                beta1[0,0,1] = beta1[0,1,0] = beta1[1,0,0] = beta_inp['xxy']
                beta1[1,0,1] = beta1[1,1,0] = beta1[0,1,1] = beta_inp['yxy']
                beta1[0,0,2] = beta1[0,2,0] = beta1[2,0,0] = beta_inp['xxz']
                beta1[0,1,2] = beta1[0,2,1] = beta1[1,0,2] = beta1[1,2,0] = beta1[2,1,0] = beta1[2,0,1] = beta_inp['yxz']
                beta1[1,1,2] = beta1[1,2,1] = beta1[2,1,1] = beta_inp['yyz']
                beta1[2,0,2] = beta1[0,2,2] = beta1[2,2,0] = beta_inp['zxz']
                beta1[2,1,2] = beta1[2,2,1] = beta1[1,2,2] = beta_inp['zyz']
                beta2 = -np.einsum('ijk,il,jm,kn->lmn', beta1, R, R, R)
                print('beta2 = \n', beta2)


                beta_vec = np.array([beta_inp['x'],beta_inp['y'],beta_inp['z']])
                res.loc[int(isomer), col_name] = np.linalg.norm(beta_vec)
res.to_csv(result_file)