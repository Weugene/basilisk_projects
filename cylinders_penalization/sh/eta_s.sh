#!/bin/bash
set -x;
#qcc -Wall -O2 cylinders_penalization.c -lm -L$BASILISK/gl  $OPENGLIBS;
for eta_s in 1e-6 1e-5 1e-4 1e-3 1e-2
do	
	echo "./a.out $eta_s >out_$eta_s 2> log_$eta_s";
	./a.out $eta_s >out_$eta_s 2> log_$eta_s;
done
