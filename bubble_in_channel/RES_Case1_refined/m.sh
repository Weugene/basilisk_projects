#!/bin/bash
set -x
echo "Usage: $0 level args"
solver=./a.out
level=$1
zlist=(0.3 1.48 14.8 2.96 0.74 1.85 3.7 9.3 0.93 1.85 3.7 7.4)
rlist=(0.075 0.075 0.075 0.075 0.15 0.15 0.15 0.15 0.3 0.3 0.3 0.3)
for i in {0..11}; do
	dir="res_${rlist[i]}_${zlist[i]}_{$level}"
	echo $dir
    (
        grep "^${level}" $dir/log_${rlist[i]}_${zlist[i]}_${level} > m_${rlist[i]}_${zlist[i]}_${level}.txt 
	)
done
find . -iname "m_*" -size  0 -print -delete
set +x
