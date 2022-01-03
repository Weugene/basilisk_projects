#!/bin/bash
set -x
set -o nounset
echo "Usage: qsub or sbatch?"
#parameters for solver
solver_parallel=./a.out
script=./run.sh

mode=$1
if [[ "$mode" == "qsub" || "$mode" == "sbatch" ]]; then
  echo "You got the right input: $mode"
else
  echo "No way. Use qsub or sbatch"
  exit
fi

Tin=(300 320 350 370 400)
Tcyl=(300)
maxlevel=(10)
iter_fp=(0 0 0 0 0 0 0 0 0 0 0 0)
ratio_dist_x=(1.1)
ratio_dist_y=(2)
ratio_front_x=-6
cyl_x=-5
Ncx=5
Ncy=8
dev_r=(0)
dev_x=(0)
dev_y=(0)
queue=cpu
Nnode=1
np_parallel=15
time_parallel=43:30:00
pmem_parallel=10GB
my_mail=evgenii.sharaborin@skoltech.ru

export PBS_O_WORKDIR=$PWD
echo "PBS_O_WORKDIR=${PBS_O_WORKDIR}"
export JLAB_SLURM_O_WORKDIR=$PWD
echo "JLAB_SLURM_O_WORKDIR=${JLAB_SLURM_O_WORKDIR}"

#	qcc $(CFLAGS) poly_heat_front.c $(CINCL) $(CLIBS) && ./a.out $(args) > out 2> log
CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  poly_heat_front.c -lm

for tc in "${Tcyl[@]}"; do
  for ti in "${Tin[@]}"; do
    for maxl in "${maxlevel[@]}"; do
      for rdx in "${ratio_dist_x[@]}"; do
        for rdy in "${ratio_dist_y[@]}"; do
          for dr in "${dev_r[@]}"; do
            for dx in "${dev_x[@]}"; do
              for dy in "${dev_y[@]}"; do
                echo "./a.out arg= ${tc} ${ti} ${maxl} ${dr} ${dx} ${dy}"
                tmp="Tin=${ti}_rdx=${rdx}"
#                tmp="Tcyl=${tc}_Tin=${ti}_maxl=${maxl}_rdx=${rdx}_rdy=${rdy}_rfx=${ratio_front_x}_Ncx=${Ncx}_Ncy=${Ncy}_dr=${dr}_dx=${dx}_dy=${dy}"
                name=${tmp}
                if [[ "$mode" == "qsub" ]]; then
                  echo "qsub"
                  string=$(qstat -a | grep -e "$name" | grep -e " R " -e " P " -e " Q ")
                  if [[ "$string" == *"$name"* ]]; then
                    echo "The submission in a queue"
                  else
                    echo "$name is not there!"
                    qsub -N "heat_${tmp}" -l nodes=${Nnode}:ppn=${np_parallel},\
                      pmem=${pmem_parallel},walltime=${time_parallel} \
                      -M ${my_mail} -m ae \
                      -F "${tc} ${ti} ${maxl} ${rdx} ${rdy} ${ratio_front_x} ${cyl_x} ${Ncx} ${Ncy} ${dr} ${dx} ${dy} ${mode}" \
                      -V ${script} &
                  fi
                else
                  echo "sbatch"
                  string=$(squeue | grep -e "$name" | grep -e " R " -e " P " -e " Q ")
                  if [[ "$string" == *"$name"* ]]; then
                    echo "The submission in a queue"
                  else
                    echo "$name is not there!"
                    sbatch --partition="${queue}" --export=ALL --job-name="${tmp}" --time=${time_parallel} \
                         --mail-user="${my_mail}" --mail-type=END,FAIL \
                         --ntasks=${np_parallel} --nodes="${Nnode}" \
                         --mem-per-cpu="${pmem_parallel}" --ntasks-per-node="${np_parallel}" \
                         --output="${tmp}-%x.%j.out" --error="${tmp}-%x.%j.log" \
                         --wrap="${script} ${tc} ${ti} ${maxl} ${rdx} ${rdy} ${ratio_front_x} ${cyl_x} ${Ncx} ${Ncy} ${dr} ${dx} ${dy} ${mode}"
                  fi
                fi
              done
            done
          done
        done
      done
    done
  done
done