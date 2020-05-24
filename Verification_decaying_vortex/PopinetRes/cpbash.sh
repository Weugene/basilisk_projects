#!/bin/bash
#echo "Usage: No arguments"
levels=(6 7 8 9 10 )

for a in "${levels[@]}"; do
        tmp=${a}_${eta}_${dtlim}
       # echo "${tmp}"
       # cat diff_${tmp}
       echo -ne $a " "
       cat res_l_${a}/log_${a} | grep "^i=" | tail -n 1
done
