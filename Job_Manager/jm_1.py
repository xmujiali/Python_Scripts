import sys, os, time, re
import subprocess as sp


'''
using this as following

    cd working_directory
    jm_1.py > out.txt 2> err.txt &

'''

if __name__ == '__main__':
    # list for storing all running processes
    procs = []
    # if such file is deleted form working_dir, whole program will stopped
    stop_code = 'delete_me_if_you_want_to_stop'
    # waiting time, it will be adjusted at runtime
    wait_sec = 10.0
    wait_sec_min, wait_sec_max = 10.0, 600.0
    alpha1,alpha2 = 7/8, 1/8

    GAUSSIAN = 'g16'
    Max_threads = 24
    free_threads = Max_threads

    working_dir = os.getcwd()
    reg = r'%nproc*?=(\d*)\s*'


    with open(stop_code, 'w') as f:
        pass

    
    while os.path.exists(working_dir + os.sep + stop_code) :
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
                    
                    succ = False
                    while not succ:
                        # if there is enough free threads
                        # submit the jobs
                        if free_threads >= nproc:
                            succ = True
                            free_threads -= nproc
                            cmd = f'{GAUSSIAN} {gjf_name}'
                            p = sp.Popen(cmd, shell=True)
                            print(cmd)
                            print('free_threads reduced to:', free_threads)
                            procs.append((p,nproc))
                            print('# procs = ', procs)
                            wait_sec = wait_sec*alpha1 + wait_sec_min*alpha2
                            print('wait_sec = ', wait_sec)
                            sys.stdout.flush()
                            time.sleep(1)
                            
                        else:
                            # if there is not enough free threads
                            # see if some previous jobs are done.
                            for item in procs:
                                p,n = item
                                if p.poll() is not None:
                                    free_threads += n
                                    procs.remove(item)
                                    print('free_threads restored to:', free_threads)
                                    print('# procs = ', procs)
                            wait_sec = wait_sec*alpha1 + wait_sec_max*alpha2
                            print('wait_sec = ', wait_sec)
                            sys.stdout.flush()
                            time.sleep(wait_sec)

                os.chdir(working_dir)
            
            
