##!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import os
import shutil
import parms
import utilities.utilities as util
from utilities.filters import filters

# goto main_path
os.chdir(parms.main_path)

wd = util.menu('Choose one working directory name:', parms.directories)
util.mkdir_safe(wd)
os.chdir(wd)

lis_templates = util.list_files(parms.main_path + os.sep + parms.gjf_template_path)
print(lis_templates)
sample_gjf_file = util.menu('Choose a gjf template:', lis_templates, custom=False)
with open(sample_gjf_file, 'r') as f:
    template = f.read()
with open('sample.gjf', 'w') as f:
    f.write(template)

ps = 'Choose a filter from following:\n' + ''.join([f"{k}\t{v['description']}\n" for k,v in filters.items()])
lis = list(filters.keys())
key = util.menu(ps, lis, custom=False)
isomers = filters[key]['filter']()


with open('isomers', 'w') as f_isomer, open('run_opt.sh', 'w') as f_run:
    for isomer in isomers:
        if isomer % parms.report_frequency == 1:
            print(f'handling isomer {isomer}.')
        f_isomer.write(str(isomer) + '\n')
        f_run.write(f'{parms.GAUSSIAN} C_{parms.total_carbons}_{isomer}.gjf\n')
        util.write_gjf(isomer, template)
        

print(f'''Job done!
Run following command lines

  cd {wd}
  bash run_opt.sh

to optimize generated isomers''')