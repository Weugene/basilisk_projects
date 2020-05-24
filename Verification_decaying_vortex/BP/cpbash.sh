#!/bin/bash
#echo "Usage: No arguments"
levels=(7 8 9 10 )
eta_s=(1e-3 1e-4 1e-5 1e-6)
dtlimiter=(0 1)

for dtlim in "${dtlimiter[@]}"; do
#echo "-------------------dtlim=${dtlim}------------"
for eta in "${eta_s[@]}"; do
#echo "================eta=${eta}==================="
for a in "${levels[@]}"; do
        tmp=${a}_${eta}_${dtlim}
       # echo "${tmp}"
       # cat diff_${tmp}
done
done
done

#echo "=============================================================="

for dtlim in "${dtlimiter[@]}"; do
#echo "-------------------dtlim=${dtlim}------------"
for a in "${levels[@]}"; do
#echo "================level=${a}==================="
for eta in "${eta_s[@]}"; do
        tmp=${a}_${eta}_${dtlim}
#        echo "${tmp}"
        echo -ne ${a} " " ${eta} " "
        cat diff_${tmp}
done
done
done

#echo "=============================================================="

for a in "${levels[@]}"; do
#echo "================level=${a}==================="
for eta in "${eta_s[@]}"; do
#echo "-------------------eta=${eta}------------"
for dtlim in "${dtlimiter[@]}"; do
        tmp=${a}_${eta}_${dtlim}
#        echo "${tmp}"
#        cat diff_${tmp}
done
done
done
