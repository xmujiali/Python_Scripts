##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import re
import os
import parms
import utilities.utilities as util

'''
假定当前目录下有failed_jobs文件,
其中记录需要重新生成结构的异构体编号,
且..\sameles\semi_opt.gjf,
是Gaussian输入文件头的模板
'''

adjacent_tab = np.zeros((parms.total_carbons+1,3), dtype=np.int32)
coords = None
with open('sample.gjf', 'r') as f_sample:
    template = f_sample.read()



# ensure_distance 和 ensure_convex
def ensure_distance(center, atom):
    global coords
    succ = True
    
    # 先计算从center指向atom的矢量(归一化)
    # 以及atom与center的距离
    vector = coords[atom] - coords[center]
    dist = np.linalg.norm(vector)
    vector /= dist

    # 处理atom是center的邻居的情形
    if atom in adjacent_tab[center]:
        if dist < parms.R_bond_min:
            coords[atom] += parms.Factor_change*vector
            coords[center] -= parms.Factor_change*vector 
            succ = False
        elif dist > parms.R_bond_max:
            coords[atom] -= parms.Factor_change*vector 
            coords[center] += parms.Factor_change*vector
            succ = False
    # 处理atom不是center邻居的情形
    else:
        if dist < parms.R_nonbond_min:
            coords[atom] += parms.Factor_change*vector
            coords[center] -= parms.Factor_change*vector 
            succ = False
    
    return succ


def ensure_convex(id_O):
    '''
    adj_tab: adjacent table
    coords: coordinates of all carbons
    id_O: currently handling carbon
    retrun : True if no modification on carbon id_0; False otherwise
    '''
    global coords
    id_A = adjacent_tab[id_O][0]
    id_B = adjacent_tab[id_O][1]
    id_C = adjacent_tab[id_O][2]

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

def clean():
    global coords
    for count in range(1,parms.Max_Cycle+1):
        succ = True
        # 取值范围 1<= center < N, i < atom <= N 
        for center in range(1,parms.total_carbons):
            for atom in range(center+1, parms.total_carbons+1):
                # 这里原来的代码移入ensure_distance
                succ = ensure_distance(center, atom) and succ
        for center in range(1,parms.total_carbons):
            # 增加 ensure_convex的步骤
            succ = ensure_convex(center) and succ
        if succ:
            # 把重心放到(0,0,0)
            coords -= coords.mean(axis=0)
            return True
    return False

# def write_gjf(isomer, coords, gjf_filename):
#     with open(gjf_filename, 'w') as f:
#         # 把头部信息写入gjf文件中
#         new_head = re.sub("\(isomer\)", str(isomer), template, re.S)
#         f.write(new_head)
#         # 把坐标写入文件
#         for c in range(1,N+1):
#             f.write(f'C\t{coords[c][0]}\t{coords[c][1]}\t{coords[c][2]}\n')
#         f.write('\n')


def reconstruct(isomer):
    global coords
    util.read_adj_tab(isomer,adjacent_tab)
    for i in range(1, parms.Max_Tries+1):
        print(f'\tTry #{i} for Isomer #{isomer}')
        coords = np.random.rand(parms.total_carbons+1,3)*8 - 3
        if clean():
            # 如果成功了，写入gjf文件并返回
            util.write_gjf(isomer, template, coords)
            return True
        # else:
            # 如果失败了，任然写入gjf文件
            # write_gjf(isomer, coords, f'C_{parms.total_carbons}_{isomer}_try{i}.gjf')
    return False

if __name__ == '__main__':
    pass
