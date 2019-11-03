#!/bin/bash
echo "Usage: ./rename.sh first step end Nstep proc"
for m in $(seq $1 $2 $3); 
do printf -v n "%04g" $((m)); 
printf -v ans "%04g" $((m+$4)); 
echo "${n} ans=${ans}"; 
for k in $(seq -f "%03g" 0 1 $5); do
	mv hrhs_${n}_n${k}.vtu hrhs_${ans}_n${k}.vtu;
done;  
mv hrhs_${n}.pvtu hrhs_${ans}.pvtu; 
done
