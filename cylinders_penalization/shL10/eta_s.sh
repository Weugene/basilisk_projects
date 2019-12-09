#!/bin/bash
set -x;
level=10
plist=(1e-6 1e-5 1e-4 1e-3 1e-2)
for eta_s in ${plist[@]};
do      
        mkdir -p res$eta_s;
        cp a.out res$eta_s/;
done
for eta_s in ${plist[@]};
do	
	echo "eta_s=$eta_s level=$level";
	./res$eta_s/a.out $eta_s $level >./res$eta_s/out_$eta_s 2> ./res$eta_s/log_$eta_s &
done

#for eta_s in ${plist[@]};
#do      
#       rm ./res$eta_s/a.out;
#done

