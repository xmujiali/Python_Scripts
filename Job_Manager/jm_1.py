import os,time,re
import multiprocessing as mp

GAUSSIAN = 'g16'
Max_threads = 24

free_threads = Max_threads
working_dir = os.getcwd()

reg = r'%nproc*?=(\d*)\s*'

def run_gauss(root, gjf_name, proc):
    os.chdir(root)
    print('\t', os.getpid(), ':',os.getcwd(), gjf_name)
    os.system(f'{GAUSSIAN} {gjf_name}')
    return proc


def done_gauss(n):
    global free_threads
    free_threads += n

if __name__ == '__main__':
    pool = mp.Pool(processes=Max_threads)
    os.chdir(working_dir)
    while True:
        for root, dirs, files in os.walk(".",):
            for gjf_name in files:
                base_name = gjf_name[:-4]
                extension = gjf_name[-4:]
                log_name = base_name + '.log'
                
                os.chdir(root)
                if '.gjf'==extension and not os.path.exists(log_name):
                    with open(gjf_name) as f:
                        nproc = 1
                        for line in f.readlines():
                            tmp = re.match(reg, line)
                            if tmp:
                                nproc = int(tmp.group(1))
                    
                    if free_threads >= nproc:
                        free_threads -= nproc
                        print(os.getpid(), ':','free_threads = ', free_threads)
                        pool.apply_async(func=run_gauss, args=(root, gjf_name, nproc), callback=done_gauss)
                        time.sleep(1)
                os.chdir(working_dir)
                

        time.sleep(300)
            
            
