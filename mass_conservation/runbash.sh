#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
rad=0.05
plist=(7 8 9 10)
#tols=(1e-4 1e-6 1e-8)
tols=(1e-8)
for tol in "${tols[@]}"; do
    for a in "${plist[@]}"; do
	    echo "./a.out arg= ${rad} ${a} ${tol}"
        tmp="r_${rad}_l_${a}_t_${tol}"
	    dir="res_${tmp}"
	    mkdir "$dir" || continue
    	cp "$solver" "$dir"
    	(
    		cd "$dir" || exit
    		mpirun -np 2 $solver "$rad" "$a" "$tol"  > "out_$tmp" 2> "log_$tmp" &
    	)
    done
done
