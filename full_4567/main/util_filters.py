#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import pandas as pd


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
