##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 18:09:27 2019
@author: Dr. Jai Li
E-mail: xmujiali@163.com
No description for this code!
"""

run_file = input('Which file you want to patition?\n')
num = input('Patiton your jobs into n parts. n = ?\n')
n = int(num)
with open(run_file, 'r') as f_run:
	lines = f_run.readlines()
	number_jobs_per_file = len(lines) // n
	for i in range(n-1):
		with open('run_partition_' + str(i) + '.sh', 'w') as part:
			for line in lines[number_jobs_per_file*i:number_jobs_per_file*(i+1)]:
				part.write(line)
	
	# handle last partition
	with open('run_partition_' + str(n-1) + '.sh', 'w') as part:
		for	line in lines[number_jobs_per_file*(n-1):]:
			part.write(line)