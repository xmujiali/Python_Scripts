import numpy as np
import parms

k1 = 1
scale = 1e-5
R_eq = 1.45

with open('failed_jobs', 'r') as f:
    for line in f.readlines():
        num = line[:-1]
        coords_file = f'..\\coords\\C_{parms.total_carbons}_{num}'
        gjf_file = f'C_{parms.total_carbons}_{num}.gjf'
        with open(coords_file, 'r') as coords, open(gjf_file, 'w') as gjf:
            coords = np.zeros((parms.total_carbons+1,3), dtype=np.float64)
            adj_table = np.zeros((parms.total_carbons+1,3),dtype=np.int32)
            text = coords.readlines()
            coords_lines = text[: parms.total_carbons]
            num = 0
            for x in coords_lines:
                num += 1
                tmp = x.split()
                x = float(tmp[1])
                y = float(tmp[2])
                z = float(tmp[3])
                coords[num][0] = x
                coords[num][1] = y
                coords[num][2] = z


            
            
            adj_lines = text[parms.total_carbons+2 : parms.total_carbons*2+2]

            for x in adj_lines:
                center, neighbors = x.split(':')
                center = int(center)
                n0,n1,n2 = neighbors.split(',')
                n0 = int(n0)
                n1 = int(n1)
                n2 = int(n2)
                print(center, n0, n1, n2)
                adj_table[center][0] = n0
                adj_table[center][1] = n1
                adj_table[center][2] = n2
            
            for center in range(1, parms.total_carbons+1):
                print(center, )
                 

