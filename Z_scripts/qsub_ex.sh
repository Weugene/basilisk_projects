#!/bin/bash
#PBS -N sleep
#PBS -l nodes=1:ppn=4
#PBS -l walltime=05:00

echo "Start date: `date`.****\n"
sleep 60
for i in {1..50}
do
   echo "Welcome $i times"
done
echo "  End date: `date`.****\n"
