#!/bin/bash
levels=( 4 5 6 7 )

for a in "${levels[@]}"; do
        tmp=${a}
       # echo "${tmp}"
	echo -ne "$a 0 "
        cat diff_${tmp} | awk '{ print $8 " " $10 }'
done
