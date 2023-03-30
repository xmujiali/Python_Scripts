import os
import multiprocessing as mp

GAUSSIAN = 'g16'
N_threads = 8


working_dir = os.getcwd()

pool = mp.Pool(processes=N_threads)
os.chdir(working_dir)
for root, dirs, files in os.walk(".",):
    for gjf_name in files:
        base_name = gjf_name[:-4]
        extension = gjf_name[-4:]
        log_name = base_name + '.log'

        if '.gjf'==extension:
            os.chdir(root)
            print(os.getcwd())
            print(f'{GAUSSIAN} {gjf_name} &')
            # os.system(f'{GAUSSIAN} {gjf_name} &')
            os.chdir(working_dir)
            




