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
import pandas as pd
import parms


def get_ring(end, old_chains):
    new_chains = []
    for chain in old_chains:
        last = chain[-1]
        if end in adj_tab[last]:
            chain.append(end)
            return chain
        for neighbor in adj_tab[last]:
            if neighbor not in chain:
                _chain = chain.copy()
                _chain.append(neighbor)
                new_chains.append(_chain)
    return get_ring(end, new_chains)


def get_all_possible_rings(end, old_chains, max_ring):
    if max_ring == parms.MAX_RING_SIZE:
        new_chains = []
        for chain in old_chains:
            last = chain[-1]
            if last == end:
                new_chains.append(chain)
        return new_chains
    
    new_chains = []
    for chain in old_chains:
        # take each chain
        last = chain[-1]
        # if this chain is already closed, we will do nothing for it
        if last == end:
            new_chains.append(chain)
        # else add one vertex for that chain
        else:
            for neighbor in adj_tab[last]:
                if neighbor not in chain:
                    _chain = chain.copy()
                    _chain.append(neighbor)
                    new_chains.append(_chain)
    return get_all_possible_rings(end, new_chains, max_ring+1)


def set_carbon_type(carbon, ring_types):
    ring_types.sort()
    _type = int(''.join((list(map(str, ring_types)))))
    car2type[carbon] = _type
    type2car[_type].append(carbon)


def set_bond_type(carbon_1, carbon_2, ring_types):
    # ensure bond (carbon_1, carbon_2) only handled once
    if carbon_1 < carbon_2:
        bond = (carbon_1, carbon_2)
        ring_types.sort()
        _type = int(''.join((list(map(str, ring_types)))))
        bond2type[bond] = _type
        type2bond[_type].append(bond)


def set_ring_type(ring):
    # make sure ring is in its unique representation and immutable
    size = len(ring)
    # transfer representation of ring from list to tuple, since list is mutable
    ring.sort()
    ring = tuple(ring)
    ring2type[ring] = size
    type2ring[size].add(ring)


def skip_lines(f, n):
    for i in range(n):
        f.readline()

#################### above are utility functions ####################

df = pd.DataFrame()
w3d = 'C%d.w3d' % parms.total_carbons
topo = 'topologies.csv'
coords_path = 'coords'
if os.path.exists(coords_path):
    shutil.rmtree(coords_path)
os.mkdir(coords_path)

# allocate some arrays to store information for each isomer
# for coding easily we use index runs from 1 to N
coords = np.zeros((parms.total_carbons+1, 3))
adj_tab = np.zeros((parms.total_carbons+1, 3), dtype='int32')


# save default value of 3 dictionaries: type2car, type2bond and size2ring
type2car_default = {}
type2bond_default = {}
type2ring_default = {}
for i in range(parms.MIN_RING_SIZE, parms.MAX_RING_SIZE+1):
    type2ring_default[i] = set()
    for j in range(i, parms.MAX_RING_SIZE+1):
        type2bond_default[int(''.join((list(map(str, [i, j])))))] = []
        for k in range(j, parms.MAX_RING_SIZE+1):
            type2car_default[int(''.join((list(map(str, [i, j, k])))))] = []


with open(w3d, 'r') as f_w3d, open(topo, 'w') as f_topo:
    # skip the first line, which is
    # >>writegraph3d<<
    # in w3d file
    skip_lines(f_w3d, 1)

    for isomer in range(1, parms.total_isomers+1):
        # skip the head line(s) of each isomer,
        # if no line need to skip, just set ISOMER_HEAD_LINE = 0 
        skip_lines(f_w3d, parms.ISOMER_HEAD_LINE)

        for c0 in range(1, parms.total_carbons+1):
            line = f_w3d.readline()
            contents = line.split()
            coords[c0, 0] = float(contents[1])  # x coordinate
            coords[c0, 1] = float(contents[2])  # y coordinate
            coords[c0, 2] = float(contents[3])  # z coordinate
            adj_tab[c0, 0] = int(contents[4])   # adjcent carbon 1
            adj_tab[c0, 1] = int(contents[5])   # adjcent carbon 2
            adj_tab[c0, 2] = int(contents[6])   # adjcent carbon 3

        # skip the tail lines, which is
        # 0
        # for each isomer
        skip_lines(f_w3d, parms.ISOMER_TAIL_LINE)

        # dicts to store carbon type, bond types and ring type
        # dict {car:car_type}
        car2type = {}
        # dict {car_type:[carbon]}
        type2car = copy.deepcopy(type2car_default)
        # dict {(c1,c2):bond_type}
        bond2type = {}
        # dict {bond_type:[(c1,c2)]}
        type2bond = copy.deepcopy(type2bond_default)
        # dict {(c1,c2,...,ck):size}
        # note we required that (c1,c2,...,ck) been ordered increasingly,
        # thus every ring has its unique representation
        ring2type = {}
        # dict {size:[(c1,c2,...,ck)]}
        type2ring = copy.deepcopy(type2ring_default)

        # starting analysis
        for c0 in range(1, parms.total_carbons+1):
            c1, c2, c3 = adj_tab[c0]
            ring12 = get_ring(c2, [[c0, c1]])
            ring12_type = len(ring12)
            ring13 = get_ring(c3, [[c0, c1]])
            ring13_type = len(ring13)
            ring23 = get_ring(c3, [[c0, c2]])
            ring23_type = len(ring23)


            if ring12_type >= (ring13_type + ring23_type - 2):
                all_possible_rings = get_all_possible_rings(c2, [[c0, c1]], max_ring=2)
                forbiden_set = (set(ring13) | set(ring23)) - set([c0,c1,c2]) 
                allowed_chains = []
                for chain in all_possible_rings:
                    if not (set(chain) & forbiden_set):
                        allowed_chains.append(chain)
                ring12 = min(allowed_chains, key=len)
                ring12_type = len(ring12)
            elif ring13_type >= (ring12_type + ring23_type - 2):
                all_possible_rings = get_all_possible_rings(c3, [[c0, c1]], max_ring=2)
                forbiden_set = (set(ring12) | set(ring23)) - set([c0,c1,c3]) 
                allowed_chains = []
                for chain in all_possible_rings:
                    if not (set(chain) & forbiden_set):
                        allowed_chains.append(chain)
                ring13 = min(allowed_chains, key=len)
                ring13_type = len(ring13)
            elif ring23_type >= (ring12_type + ring13_type - 2):
                all_possible_rings = get_all_possible_rings(c3, [[c0, c2]], max_ring=2)
                forbiden_set = (set(ring12) | set(ring13)) - set([c0,c2,c3]) 
                allowed_chains = []
                for chain in all_possible_rings:
                    if not (set(chain) & forbiden_set):
                        allowed_chains.append(chain)
                ring23 = min(allowed_chains, key=len)
                ring23_type = len(ring23)


            set_carbon_type(c0, [ring12_type, ring13_type, ring23_type])
            set_bond_type(c0, c1, [ring12_type, ring13_type])
            set_bond_type(c0, c2, [ring12_type, ring23_type])
            set_bond_type(c0, c3, [ring13_type, ring23_type])
            set_ring_type(ring12)
            set_ring_type(ring13)
            set_ring_type(ring23)
        
        
        # save topological statistics
        for t, cars in type2car.items():
            df.loc[isomer, t] = len(cars)
        for t, bonds in type2bond.items():
            df.loc[isomer, t] = len(bonds)
        for t, rings in type2ring.items():
            df.loc[isomer, t] = len(rings)


        # save coordinates and other infomation
        filename = 'C_%d_%d' % (parms.total_carbons, isomer)
        with open(os.path.join(coords_path, filename), 'w') as f_coord:
            for c0 in range(1, parms.total_carbons+1):
                f_coord.write('C\t%9.5f\t%9.5f\t%9.5f\n' % (coords[c0][0], coords[c0][1], coords[c0][2]))
            f_coord.write('\nadjencent table:\n')
            for c0 in range(1, parms.total_carbons+1):
                f_coord.write('%d:%d,%d,%d\n' % (c0, adj_tab[c0][0], adj_tab[c0][1], adj_tab[c0][2]))
            f_coord.write('\ncar2type:\n')
            f_coord.write(str(car2type) + '\n')
            f_coord.write('\ntype2car:\n')
            f_coord.write(str(type2car) + '\n')
            f_coord.write('\nbond2type:\n')
            f_coord.write(str(bond2type) + '\n')
            f_coord.write('\ntype2bond:\n')
            f_coord.write(str(type2bond) + '\n')
            f_coord.write('\nring2type:\n')
            f_coord.write(str(ring2type) + '\n')
            f_coord.write('\ntype2ring:\n')
            f_coord.write(str(type2ring) + '\n')

df.to_csv(topo)
