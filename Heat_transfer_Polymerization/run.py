#!/usr/bin/env python3

import yaml
import sys
import subprocess
import numpy as np
import os

file = sys.argv[1]
myenv = os.environ.copy()
print('LD_LIBRARY_PATH:', myenv['LD_LIBRARY_PATH'], ' C_INCLUDE_PATH:', myenv['C_INCLUDE_PATH'])
Uin = np.array([0.005, 0.1,  0.2])
Tin = np.array([300, 350])
T_solid = np.array([300, 350, 400])
#Uin = [0.01]
#Tin = [300]
#T_solid = [350]

parallel = 10
Nnode = 1
queue = "cpu"
time_parallel = "50:00:00"
pmem_parallel = "10GB"
my_mail = "evgenii.sharaborin@skoltech.ru"
outfile = 'out'
logfile = 'log'
maxTime = 10000000
command = 'rm a.out; make poly_with_yaml'
print(command + '...')
ret = subprocess.run(command, capture_output=True, shell=True)
print("finished")
for tin in Tin:
    for t_solid in T_solid:
        for uin in Uin:
            dirn = f'Tin={tin}_Tsolid={t_solid}_Uin={uin}'
            command = f'mkdir {dirn};'
            command += f'cd {dirn};'
            command += f'cp ../{file} . ;'
            command += f'cp ../a.out . ;'
            subprocess.run(command, capture_output=True, shell=True)
            print(command)
            with open(os.path.join(dirn, file)) as f:
                data = yaml.safe_load(f)
            data['name'] = str(dirn)
            data['dimensional_vars']['Uin'] = float(uin)
            data['dimensional_vars']['Tin'] = float(tin)
            data['dimensional_vars']['Tam'] = float(tin)
            data['dimensional_vars']['T_solid'] = float(t_solid)

            with open(os.path.join(dirn, file), 'w') as outf:
                yaml.dump(data, outf, default_flow_style=False)


            # command = f'cd {dirn} && mpirun -np {parallel} ./a.out {file} > {outfile} 2> {logfile} &'

            command = f'''cd {dirn} && sbatch --partition="{queue}" \
            --export=ALL --job-name="{dirn}" \
            --time={time_parallel} \
            --mail-user="{my_mail}" --mail-type=END,FAIL \
            --ntasks={parallel} --nodes="{Nnode}" \
            --mem-per-cpu="{pmem_parallel}" \
            --output="out_{dirn}" --error="log_{dirn}" \
            --wrap="../sbatch {dirn} {file}"'''
            print(command + '...')
            subprocess.run(command, shell=True)
            print("finished")

subprocess.run(f"sleep {maxTime}; echo 2", capture_output=True, shell=True)
