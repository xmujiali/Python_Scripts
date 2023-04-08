import re
if __name__ == '__main__':
    N = 30
    failed_jobs = 'failed_jobs'
    reg_exp = '#'
    routine_section = input('Input the new routine section:\n') + '\n'

    
    with open(failed_jobs,'r') as f:
        for line in f.readlines():
            # 得到异构体的编号
            isomer = int(line.strip())
            print(f'Handling Isomer: {isomer}')
            gjf_name = f'C_{N}_{isomer}.gjf'
            with open(gjf_name, 'r') as gjf:
                lines = gjf.readlines()
                
            for i,x in enumerate(lines):
                if re.match(reg_exp, x):
                     lines[i] = routine_section
                     
            with open(gjf_name, 'w') as gjf:
                 gjf.writelines(lines)
                 gjf.write('\n')
                



