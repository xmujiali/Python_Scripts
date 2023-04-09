##!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import copy
import shutil
import numpy as np
import numpy as np
import pandas as pd
import configparser

# initialize paramters
parms = configparser.ConfigParser()
parms.read('.\setting.ini', encoding='UTF-8')

MIN_RING_SIZE = parms.getint('topo', 'MIN_RING_SIZE')
MAX_RING_SIZE = parms.getint('topo', 'MAX_RING_SIZE')
total_carbons = parms.getint('topo', 'total_carbons')
total_isomers = parms.getint('topo', 'total_isomers')
ISOMER_HEAD_LINE = parms.getint('w3d', 'ISOMER_HEAD_LINE')
ISOMER_TAIL_LINE = parms.getint('w3d', 'ISOMER_TAIL_LINE')

def get_ring(v0,v1,v2,v3):
    '''
    we need to distinguish ring and cycle:
    ring is defind as polygon in fullerene surface;
    while cycle, as a graph theoretical concept, is defined as simple path with initial vertex and final vertex are identical.
    ring is simply one cycle. But a cycle need not to be a ring, which could be boundary of a multi-ring fused area.

    To make sure a cycle is a ring, it must satisfies two conditions:
    (A) (NO EDGE) For any vertex, only 2 neighbor in the ring, NOT 3.
    otherwise， there will be vertices or edge inside the cycle, which means it is NOT a ring.
    (B) (NO VERTEX) All the rest vertices, except ones in the ring, should be connected. Namely, from any one vertex, we collect its neighbors and neighbors of neighbors...and so on (except the vertices in the cycle), we will finally get the all the vertices except the ones in the ring.
        
    '''
    # result stores cycles found by get_cycles
    # tmp_chains stores uncompleted chains, which is not cycles, but may be later.
    cycles = {}
    for i in range(MIN_RING_SIZE,MAX_RING_SIZE+1):
        cycles[i] = []
    tmp_chains = [[v0,v1]]
    end = v2
    # v3 is 3rd neighbor of v0, which is of course not in the ring. 
    never_include = v3

    def get_cycles(current_size):
        '''
        get cycles will generate all cycles like [v0,v1,...,vN,v0], 
        with N < MAX_RING_SIZE
        '''
        nonlocal tmp_chains
        new_chains = []
        for chain in tmp_chains:
            last = chain[-1]
            # take each chain, add one one neighbor of last vertex for this chain
            for neighbor in adj_tab[last]:
                if neighbor not in chain:
                    _chain = chain.copy()
                    _chain.append(neighbor)
                    # _chain is already a cycle, we will store it to cycles
                    if neighbor == end:
                        _chain.sort()
                        cycles[len(_chain)].append(_chain)
                    # otherwise we will store it to new_chains
                    else:
                        new_chains.append(_chain)
        
        current_size += 1
        if current_size < MAX_RING_SIZE:
            # update tmp_chains and call get_cycles() recursively
            tmp_chains = new_chains
            get_cycles(current_size)
        
        # return explicitly
        return

    def cond_A(cycle):
        for v in cycle:
            if len(set(adj_tab[v]) & set(cycle)) == 3:
                return False
        return True

    def cond_B(cycle):
        cycle = set(cycle)
        explored = set()
        new_ones = set([never_include])
        while len(new_ones):
            new_new_ones = set()
            for x in new_ones:
                new_new_ones = new_new_ones | ( set(adj_tab[x]) - cycle - explored)
            explored = explored | new_ones
            new_ones = new_new_ones
        if len(explored) + len(cycle) == total_carbons:
            return True
        else:
            return False

    get_cycles(2)
    result = []
    for xs in cycles.values():
        for x in xs: 
            if cond_A(x) and cond_B(x):
                x.sort()
                result.append(x)
    if len(result) != 1:
        print('#isomer = ', isomer)
        print('there is something wrong with parmater get_ring(%d, %d, %d, %d)' % (v0, v1, v2, v3))
    return result[0]

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

def count_kekule(sites):
    length = len(sites)
    if length == 0:
        return 1
    
    num_kekule = 0
    last = sites[-1]
    # i在sites且是last的邻居
    nbrs = set(sites) & set(adj_tab[last])
    for i in nbrs:
        lis_ = sites.copy()[:-1]
        lis_.remove(i)
        num_kekule += count_kekule(lis_)
    return num_kekule

def skip_lines(f, n):
    for i in range(n):
        f.readline()

#################### above are utility functions ####################

df_topo = pd.DataFrame()
df_Kekule = pd.DataFrame()
df_Huckel = pd.DataFrame()
w3d = f'C{total_carbons}.w3d' 
topo = 'topologies.csv'
kekule = 'kekule.csv'
huckel = 'huckel.csv'
coords_path = 'coords'
if os.path.exists(coords_path):
    shutil.rmtree(coords_path)
os.mkdir(coords_path)

# allocate arrays to store information for each isomer
# for coding easily, we use index runs from 1 to N
coords = np.zeros((total_carbons+1, 3))
adj_tab = np.zeros((total_carbons+1, 3), dtype='int32')
adj_matrix = np.zeros((total_carbons, total_carbons), dtype='int32')
all_carbons = list(range(1, total_carbons+1))


# save default value of 3 dictionaries: type2car, type2bond and size2ring
type2car_default = {}
type2bond_default = {}
type2ring_default = {}
for i in range(MIN_RING_SIZE, MAX_RING_SIZE+1):
    type2ring_default[i] = set()
    for j in range(i, MAX_RING_SIZE+1):
        type2bond_default[int(''.join((list(map(str, [i, j])))))] = []
        for k in range(j, MAX_RING_SIZE+1):
            type2car_default[int(''.join((list(map(str, [i, j, k])))))] = []


with open(w3d, 'r') as f_w3d:
    # skip the first line, which is
    # >>writegraph3d<<
    # in w3d file
    skip_lines(f_w3d, 1)

    for isomer in range(1, total_isomers+1):
        # skip the head line(s) of each isomer,
        # if no line need to skip, just set ISOMER_HEAD_LINE = 0 
        skip_lines(f_w3d, ISOMER_HEAD_LINE)
        adj_matrix[:, :] = 0 

        for c0 in range(1, total_carbons+1):
            line = f_w3d.readline()
            contents = line.split()
            coords[c0, 0] = float(contents[1])  # x coordinate
            coords[c0, 1] = float(contents[2])  # y coordinate
            coords[c0, 2] = float(contents[3])  # z coordinate

            nb1, nb2, nb3 = int(contents[4]), int(contents[5]), int(contents[6])
            adj_tab[c0, 0] = nb1   # adjcent carbon 1
            adj_tab[c0, 1] = nb2   # adjcent carbon 2
            adj_tab[c0, 2] = nb3   # adjcent carbon 3
            # note the index of adj_matrix runs from 0 to N-1
            adj_matrix[c0-1, nb1-1] = 1
            adj_matrix[c0-1, nb2-1] = 1
            adj_matrix[c0-1, nb3-1] = 1

        # skip the tail lines, which is
        # 0
        # for each isomer
        skip_lines(f_w3d, ISOMER_TAIL_LINE)

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
        for c0 in range(1, total_carbons+1):
            c1, c2, c3 = adj_tab[c0]
            ring12 = get_ring(c0, c1, c2, c3)
            ring12_type = len(ring12)
            ring13 = get_ring(c0, c1, c3, c2)
            ring13_type = len(ring13)
            ring23 = get_ring(c0, c2, c3, c1)
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
            df_topo.loc[isomer, t] = len(cars)
        for t, bonds in type2bond.items():
            df_topo.loc[isomer, t] = len(bonds)
        for t, rings in type2ring.items():
            df_topo.loc[isomer, t] = len(rings)


        # save coordinates and other infomation
        filename = 'C_%d_%d' % (total_carbons, isomer)
        with open(os.path.join(coords_path, filename), 'w') as f_coord:
            for c0 in range(1, total_carbons+1):
                f_coord.write('C\t%9.5f\t%9.5f\t%9.5f\n' % (coords[c0][0], coords[c0][1], coords[c0][2]))
            f_coord.write('\nadjencent table:\n')
            for c0 in range(1, total_carbons+1):
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
    
        # Huckel
        evs = np.linalg.eigvalsh(-adj_matrix)
        evs.sort()
        for i in range(len(evs)):
            orb = 'orb' + str(i+1)
            df_Huckel.loc[isomer, orb] = evs[i]

        # Kekule
        k = count_kekule(all_carbons)
        df_Kekule.loc[isomer, 'K'] = k
        df_Kekule.loc[isomer, 'ln_K'] = np.log(k)

df_topo.to_csv(topo)
df_Huckel.to_csv(huckel)
df_Kekule.to_csv(kekule)
