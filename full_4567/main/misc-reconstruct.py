#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import re


'''
本程序要在C30\semi_opt\目录下运行！
'''


# 一些参数设定
sample_gjf_file = '..\\samples\\semi_opt.gjf'
failed_jobs = 'failed_jobs'
Max_Tries = 12
Max_Cycle = 10000
R0 = 2.2
R_max = 1.6
R_min = 1.4
Factor_change = 0.02

N = 30
# 参数结束




# 增加两个新的函数
# ensure_distance 和 ensure_convex
def ensure_distance(adj_tab, coords, center, atom):
    succ = True

    # 先计算从center指向atom的矢量(归一化)
    # 以及atom与center的距离
    vector = coords[atom] - coords[center]
    dist = np.linalg.norm(vector)
    vector /= dist

    # 处理atom是center的邻居的情形
    if atom in adj_tab[center]:
        if dist < R_min:
            coords[atom] += Factor_change*vector
            coords[center] -= Factor_change*vector 
            succ = False
        elif dist > R_max:
            coords[atom] -= Factor_change*vector 
            coords[center] += Factor_change*vector
            succ = False
    # 处理atom不是center邻居的情形
    else:
        if dist < R0:
            coords[atom] += Factor_change*vector
            coords[center] -= Factor_change*vector 
            succ = False
    
    return succ


def ensure_convex(adj_tab, coords, id_O):
    '''
    参数说明如下：
    adj_tab: 邻接表
    coords: 一组包含所有原子的坐标，
    id_O, 待处理原子O的编号
    返回值，没有更改id_O原子的坐标返回True，否则返回False
    '''
    id_A = adj_tab[id_O][0]
    id_B = adj_tab[id_O][1]
    id_C = adj_tab[id_O][2]

    # vec1 = point_C - point_A
    vec1 = coords[id_C] - coords[id_A]

    # vec2 = point_B - point_A
    vec2 = coords[id_B] - coords[id_A]

    # vec = point_O - point_A
    vec = coords[id_O] - coords[id_A]

    # vec_U = bary_center - point_A
    bary_center = coords.mean(axis=0)
    vec_U = bary_center - coords[id_A]

    w = np.cross(vec1, vec2)
    w /= np.linalg.norm(w)
    vec_orth =  np.dot(vec, w) * w

    # 条件成立说明vec_U, vec 在垂直平面方向的分量同号，
    # 即O和barycenter在平面的同一侧
    if np.dot(vec_U,w) * np.dot(vec, w) > 0:
        vec_ =  vec - 2 * vec_orth

        # point_O = point_A + vec_ 
        coords[id_O] = coords[id_A] + vec_ 
        return False

    return True

def clean(adj_tab, coords):
    count = 1
    while count <= Max_Cycle:
        succ = True
        # 取值范围 1<= center < N, i < atom <= N 
        for center in range(1,N):
            for atom in range(center+1, N+1):
                # 这里原来的代码移入ensure_distance
                succ = ensure_distance(adj_tab, coords, center, atom) and succ
            
        for center in range(1,N):
            # 增加 ensure_convex的步骤
            succ = ensure_convex(adj_tab, coords, center) and succ

        if succ:
            # 把重心放到(0,0,0)
            print(coords.mean(axis=0))
            coords -= coords.mean(axis=0)
            return True
        count += 1
        # print(f'\t\tcount = {count}')
    
    return False

def write_gjf(iso, coords, gjf_filename):
    with open(gjf_filename, 'w') as f:
        # 把头部信息写入gjf文件中
        new_head = re.sub("\(isomer\)", str(iso), head, re.S)
        f.write(new_head)
        # 把坐标写入文件
        for c in range(1,N+1):
            f.write(f'C\t{coords[c][0]}\t{coords[c][1]}\t{coords[c][2]}\n')
        f.write('\n')

def reconstruct(iso, adj_tab):
    for i in range(1, Max_Tries+1):
        print(f'\tTry #{i} for Isomer #{isomer}')
        coords = np.random.rand(N+1,3)*6 - 3
        if clean(adj_tab, coords):
            # 如果成功了，写入gjf文件并返回
            write_gjf(iso, coords, f'C_{N}_{iso}.gjf')
            return True
        else:
            # 如果失败了，任然写入gjf文件
            write_gjf(iso, coords, f'C_{N}_{iso}_try{i}.gjf')
    return False




with open(sample_gjf_file, 'r') as f_sample:
    head = f_sample.read()
    
with open(failed_jobs,'r') as f:
    for line in f.readlines():
        # 得到异构体的编号
        isomer = int(line.strip())
        print(f'Handling Isomer: {isomer}')
        with open(f'..\\coords\\C_{N}_{isomer}', 'r') as g:
            adjacent_tab = np.zeros((N+1,3), dtype=np.int32)
            
            for line in g.readlines()[N+2:2*N+2]:
                center, neighbors = line.split(':')
                center = int(center)
                nbr1,nbr2,nbr3 = neighbors.split(',')
                nbr1 = int(nbr1)
                nbr2 = int(nbr2)
                nbr3 = int(nbr3)
                adjacent_tab[center] = np.array([nbr1,nbr2,nbr3])
            # print('adj_tab =\n', adjacent_tab)
            
            if reconstruct(isomer, adjacent_tab):
                print(f'reconstruction of Isomer #{isomer} is success!')
            else:
                print(f'reconstruction of Isomer #{isomer} is failed!')




