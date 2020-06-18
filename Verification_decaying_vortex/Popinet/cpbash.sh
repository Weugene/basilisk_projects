#!/bin/bash
echo "Usage: No arguments"
levels=( 4 5 6 7 )
for a in "${levels[@]}"; do
	tmp=${a}
	dir="res_l_${tmp}"
	echo "${dir}"
	(
		cd "$dir" || exit
		grep "i=" "log_${tmp}" |tail -n 1 > diff_${tmp}
		cp diff_${tmp} ../
	)
done

for a in "${levels[@]}"; do
        tmp=${a}
        echo "${tmp}"
        cat diff_${tmp}
done
