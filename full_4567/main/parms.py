#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 21:37:22 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""
import os
import pandas as pd

# configuration start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

total_carbons = 30
total_isomers = 26639
# total_carbons = 60
# total_isomers = 1812
charge = 0
MIN_RING_SIZE, MAX_RING_SIZE = 4, 7
ISOMER_HEAD_LINE = 0
# ISOMER_HEAD_LINE = 1
ISOMER_TAIL_LINE = 1
max_Active_Sites = 6
GAUSSIAN = 'g16'



# following are filters
# each filter take one DataFrame containing all topological information
# and return a list of isomer numbers
def choose_all():
        file = input('Give me csv filename[Enter for ./topologies.csv]:')
        if file == '':
                file = '.' + os.sep + 'topologies.csv'
        df = pd.read_csv(file, index_col=0)
        return df.index.tolist()

def larger_than():
        file = input('Give me csv filename[Enter for ./result.csv]:')
        if file == '':
                file = '.' + os.sep + 'result.csv'
        df = pd.read_csv(file, index_col=0)
        print(df.columns)
        col_name = input('col_name = ')
        print(df[col_name].describe())
        value = float(input('Only isomers with value larger than you give will be considered.\nvalue = '))
        df = df[df[col_name] > value]
        return df.index.tolist()

def less_than():
        file = input('Give me csv filename[Enter for ./result.csv]:')
        if file == '':
                file = '.' + os.sep + 'result.csv'
        df = pd.read_csv(file, index_col=0)
        print(df.columns)
        col_name = input('col_name = ')
        print(df[col_name].describe())
        value = float(input('Only isomers with value less than you give will be considered.\nvalue = '))
        df = df[df[col_name] < value]
        return df.index.tolist()


filters = {
        'all' : {'description' : 'choose all isomers.',
                'filter' : choose_all},
        'GT' : {'description' : 'returns isomers with value larger than you give.',
                'filter' : larger_than},
        'LT' : {'description' : 'returns isomers with value less than you give.',
                'filter' : less_than},
}


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< configuration end
