##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 19:04:00 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

import os
import copy
import shutil
import numpy as np
import parms


def get_ring(target: int, old_chains: [[int]]):
    """
    :parm target: carbon serial of end point of a ring
    :parm old_chains: a list of carbon chain (a list), the returned ring will
        be extension of one of these chains.
    :return: a smallest carbon ring (a list), the carbons in the list have
        already sorted by its serial
    This function extend every chain in old_chains by one carbon at its tail
        util one of the new chain reached the target carbon.
    if user want to find a smallest ring [1,2, ... , 6], note 6 is connected
        back to 1, call this function as get_ring(6, [[1,2]]).
    """
    new_chains = []
    for chain in old_chains:
        last = chain[-1]
        next_last = chain[-2]
        if target in adj_tab[last]:
            chain.append(target)
            chain.sort()
            return tuple(chain)
        for neighbor in adj_tab[last]:
            if neighbor != next_last:
                _chain = chain.copy()
                _chain.append(neighbor)
                new_chains.append(_chain)
    return get_ring(target, new_chains)


def set_carbon_type(carbon, ring_size_lis):
    """
    TO BE DONE!!!
    """
    ring_size_lis.sort()
    _type = int(''.join((list(map(str, ring_size_lis)))))
    car2type[carbon] = _type
    type2car[_type].append(carbon)


def set_bond_type(c_a, c_b, ring_size_lis):
    """
    TO BE DONE!!!
    """
    # make sure bond (c_a, c_b) only handled once
    if c_a > c_b:
        return
    bond = (c_a, c_b)
    ring_size_lis.sort()
    _type = int(''.join((list(map(str, ring_size_lis)))))
    bond2type[bond] = _type
    type2bond[_type].append(bond)


def set_ring_size(ring):
    """
    TO BE DONE!!!
    Note: make sure ring is already in its unique representation and immutable
    """
    size = len(ring)
    # transfer representation of ring from list to tuple, since list is mutable
    ring = tuple(ring)
    ring2size[ring] = size
    size2ring[size].add(ring)


def skip_lines(f: object, n: int):
    for i in range(n):
        f.readline()

#################### above are utility functions ####################


w3d = 'C%d.w3d' % parms.total_carbons
topo = 'topologies.csv'
# multiplicity = parms.charge % 2 + 1

# allocate some arrays to store information for each isomer
coords = np.zeros((parms.total_carbons+1, 3))
adj_tab = np.zeros((parms.total_carbons+1, 3), dtype='int32')


# save default value of 3 dictionaries: type2car, type2bond and size2ring
type2car_default = {}
type2bond_default = {}
size2ring_default = {}
for i in range(parms.MIN_RING_SIZE, parms.MAX_RING_SIZE+1):
    size2ring_default[i] = set()
    for j in range(i, parms.MAX_RING_SIZE+1):
        type2bond_default[int(''.join((list(map(str, [i, j])))))] = []
        for k in range(j, parms.MAX_RING_SIZE+1):
            type2car_default[int(''.join((list(map(str, [i, j, k])))))] = []


with open(w3d, 'r') as w3d, open(topo, 'w') as topo:
    # make directory coords
    coords_path = 'coords'
    if os.path.exists(coords_path):
        shutil.rmtree(coords_path)
    os.mkdir(coords_path)

    # write the head line in topologies.csv
    head_line = '#Isomer'
    for t in type2car_default:
        head_line += ',#C%d' % t
    for t in type2bond_default:
        head_line += ',#B%d' % t
    for size in size2ring_default:
        head_line += ',#R%d' % size
    head_line += '\n'
    # head_line += '#Kekule,E_HMO\n'
    topo.write(head_line)

    # skip the first line in w3d file, like
    # >>writegraph3d<<
    skip_lines(w3d, 1)

    for isomer in range(1, parms.total_isomers+1):
        # skip the head line(s) of each isomer, like
        # nmrdata      269   C2    25 xcarbon  2 x
        skip_lines(w3d, parms.ISOMER_HEAD_LINE)

        for car in range(1, parms.total_carbons+1):
            line = w3d.readline()
            contents = line.split()
            coords[car, 0] = float(contents[1])  # x coordinate
            coords[car, 1] = float(contents[2])  # y coordinate
            coords[car, 2] = float(contents[3])  # z coordinate
            adj_tab[car, 0] = int(contents[4])   # adjcent carbon 1
            adj_tab[car, 1] = int(contents[5])   # adjcent carbon 2
            adj_tab[car, 2] = int(contents[6])   # adjcent carbon 3

        # skip the tail lines of each isomer, like
        # 0
        skip_lines(w3d, parms.ISOMER_TAIL_LINE)

        # dicts to store carbon type, bond types and ring size
        # dict {car:car_type}
        car2type = {}
        # dict {car_type:[carbon]}
        type2car = copy.deepcopy(type2car_default)
        # dict {(c1,c2):bond_type}
        bond2type = {}
        # dict {bond_type:[(c1,c2)]}
        type2bond = copy.deepcopy(type2bond_default)
        # dict {(c1,c2,...,ck):size}
        ring2size = {}
        # dict {size:[(c1,c2,...,ck)]}
        size2ring = copy.deepcopy(size2ring_default)

        # starting analsys
        for car in range(1, parms.total_carbons+1):
            c1, c2, c3 = adj_tab[car]
            ring12 = get_ring(c2, [[car, c1]])
            ring_size12 = len(ring12)
            ring13 = get_ring(c3, [[car, c1]])
            ring_size13 = len(ring13)
            ring23 = get_ring(c3, [[car, c2]])
            ring_size23 = len(ring23)

            set_carbon_type(car, [ring_size12, ring_size13, ring_size23])
            set_bond_type(car, c1, [ring_size12, ring_size13])
            set_bond_type(car, c2, [ring_size12, ring_size23])
            set_bond_type(car, c3, [ring_size13, ring_size23])
            set_ring_size(ring12)
            set_ring_size(ring13)
            set_ring_size(ring23)

        # save topological statistics
        line = '%d' % isomer
        for t, cars in type2car.items():
            line += ',%d' % len(cars)
        for t, bonds in type2bond.items():
            line += ',%d' % len(bonds)
        for size, rings in size2ring.items():
            line += ',%d' % len(rings)
        line += '\n'
        # line += ',%d,%.6f\n' % (1, 1.0)
        topo.write(line)

        # save coordinates and other infomation
        filename = 'C_%d_%d' % (parms.total_carbons, isomer)
        with open(os.path.join(coords_path, filename), 'w') as coord:
            for car in range(1, parms.total_carbons+1):
                coord.write('C\t%9.5f\t%9.5f\t%9.5f\n'
                            % (coords[car][0], coords[car][1], coords[car][2]))
            coord.write('\n****************************************:\n')
            coord.write('\nadjencent table:\n')
            for car in range(1, parms.total_carbons+1):
                coord.write('%d:%d,%d,%d\n'
                            % (car, adj_tab[car][0], adj_tab[car][1], adj_tab[car][2]))
            coord.write('\n****************************************:\n')
            coord.write('\ncar2type:\n')
            coord.write(str(car2type) + '\n')
            coord.write('\n****************************************:\n')
            coord.write('\ntype2car:\n')
            coord.write(str(type2car) + '\n')
            coord.write('\n****************************************:\n')
            coord.write('\nbond2type:\n')
            coord.write(str(bond2type) + '\n')
            coord.write('\n****************************************:\n')
            coord.write('\ntype2bond:\n')
            coord.write(str(type2bond) + '\n')
            coord.write('\n****************************************:\n')
            coord.write('\nring2size:\n')
            coord.write(str(ring2size) + '\n')
            coord.write('\n****************************************:\n')
            coord.write('\nsize2ring:\n')
            coord.write(str(size2ring) + '\n')
