#!/bin/bash
set -x
echo "Usage: $0 level args"
solver=./a.out #Rad Ca Re maxlevel Rmu
level=$1
Re=0.1
calist=(0.005 0.01 0.015 0.02 0.035 0.045 0.1 0.25)
rmulist=(0.125 0.4 1)
rad=0.45

for ca in "${calist[@]}"; do
    for rmu in "${rmulist[@]}"; do
        name="${rad}_${ca}_${Re}_${level}_${rmu}"
	    dir="res_${name}"
	    echo $dir
        (
            grep "^${level}" $dir/log_${name} > m_${name}.txt 
	    )
    done
done
find . -iname "m_*" -size  0 -print -delete
set +x
