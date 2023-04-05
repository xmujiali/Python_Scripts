##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 19:04:00 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import numpy as np
import pandas as pd
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

df = pd.read_csv('topologies.csv', index_col=0, dtype=int)
Vertices, Edges, Faces = [], [], ['4', '5', '6', '7']
for f in df.columns:
    if len(f)==3:
        # 处理顶点
        Vertices.append(f)
    elif len(f)==2:
        # 处理边
        Edges.append(f)

for i,d in df.iterrows():
    logging.debug('handling isomer # %5d：' %i)
    # Euler公式
    _F, _E, _V = 0, 0, 0
    for f in Faces:
        _F += d[f]
    for e in Edges:
        _E += d[e]
    for v in Vertices:
        _V += d[v]
    if _V - _E + _F != 2:
        logging.debug('There is a contradiction for Euler Theorem!')
        logging.debug('V,E,F = ', _V, _E, _F)

    # 面和顶点的关系
    for f in Faces:
        LHS = int(f) * d[f]
        LHS_str = str(int(f)) + ' * (F' + f + '=' + str(d[f]) + ')'
        RHS = 0
        RHS_str = ''
        for v in Vertices:
            coefficient = v.count(f)
            if coefficient > 0:
                RHS +=  coefficient * d[v]
                RHS_str += ' + ' + str(coefficient) + '* (V' + v + '=' + str(d[v]) + ')'
        if LHS != RHS:
            logging.debug('-'*20)
            logging.debug('There is a problem for F-V relation')
            logging.debug(LHS_str + ' = ' + RHS_str)
    

    # 面和边的关系
    for f in Faces:
        LHS = int(f) * d[f]
        LHS_str = str(int(f)) + ' * (F' + f + '=' + str(d[f]) + ')'
        RHS = 0
        RHS_str = ''
        for e in Edges:
            coefficient = e.count(f)
            if coefficient > 0:
                RHS +=  coefficient * d[e]
                RHS_str += ' + ' + str(coefficient) + '* (E' + e + '=' + str(d[e]) + ')'
        if LHS != RHS:
            logging.debug('-'*20)
            logging.debug('There is a problem for F-E relation')
            logging.debug(LHS_str + ' = ' + RHS_str)

    # 边和顶点的关系
    for e in Edges:
        f1, f2 = e[0], e[1]
        LHS = 2 * d[e]
        LHS_str = str(2) + ' * (E' + e + '=' + str(d[e]) + ')'
        RHS = 0
        RHS_str = ''
        for v in Vertices:
            if f1 == f2:
                tmp = v.count(f1)
                if tmp == 3:
                    coefficient = 3
                elif tmp == 2:
                    coefficient = 1
                else:
                    coefficient = 0
            else:
                coefficient = v.count(f1) * v.count(f2)
            
            if coefficient > 0:
                RHS +=  coefficient * d[v]
                RHS_str += ' + ' + str(coefficient) + '* (V' + v + '=' + str(d[v]) + ')'
        if LHS != RHS:
            logging.debug('-'*20)
            logging.debug('There is a problem for E-V relation')
            logging.debug(LHS_str + ' = ' + RHS_str)



