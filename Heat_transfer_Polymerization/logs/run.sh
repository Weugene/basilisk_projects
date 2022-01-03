#!/bin/bash
set -x
set -o nounset
echo "Usage: Tcyl Tin maxlevel dr dx dy"

Tcyl=$1
Tin=$2
maxlevel=$3
ratio_Rbmin=0.166666
ratio_Rbmax=1.
ratio_dist_x=$4
ratio_dist_y=$5
ratio_front_x=$6
cyl_x=$7
Nb=0
Ncx=$8
Ncy=$9
TOLERANCE_P=1e-6
TOLERANCE_V=1e-6
TOLERANCE_T=1e-6
Htr=355000
Arrhenius_const=80600
Ea_by_R=7697 #already divided by 8.31
dev_r=${10}
dev_x=${11}
dev_y=${12}
mode=${13}


tmp="Tcyl=${Tcyl}_Tin=${Tin}_maxl=${maxlevel}_rdx=${ratio_dist_x}_rdy=${ratio_dist_y}_rfx=${ratio_front_x}_Ncx=${Ncx}_Ncy=${Ncy}_dr=${dev_r}_dx=${dev_x}_dy=${dev_y}"
mkdir ${tmp} || continue
cd ${tmp} || continue
cp ../a.out .

#find out iter_fp
myarr=($(grep "<DataSet timestep=" heat_pol.pvd ))
NN=${#myarr[@]}
if (($NN > 0)); then
	last_iter_fp=$(echo ${myarr[${NN}-1]} | sed 's/file="res\/heat\_pol\_0\_//'| sed 's/\.pvtu"\/>//' | bc)
fi
if ((${last_iter_fp} > 0)); then
  iter_fp=${last_iter_fp}
else
  iter_fp=0
fi

args="${Tcyl} ${Tin} ${maxlevel} ${iter_fp} ${ratio_Rbmin} ${ratio_Rbmax} ${ratio_dist_x} ${ratio_dist_y} ${ratio_front_x} \
${cyl_x} ${Nb} ${Ncx} ${Ncy} ${TOLERANCE_P} ${TOLERANCE_V} ${TOLERANCE_T} \
${Htr} ${Arrhenius_const} ${Ea_by_R} ${dev_r} ${dev_x} ${dev_y}"

echo "args:${args}"

echo "BASILISK=$BASILISK"


#save the previous heat_pol.pvd and gather all in result.pvd

lines_in_prev=$(head -n -2 heat_pol_prev.pvd | wc -l)
if (( $lines_in_prev > 0 )); then
  head -n -2 heat_pol_prev.pvd > result.pvd
  tail -n +$(($lines_in_prev + 1)) heat_pol.pvd >> result.pvd
else
  cp heat_pol.pvd heat_pol_prev.pvd || continue
  cp heat_pol.pvd result.pvd || continue
fi
cp heat_pol.pvd heat_pol_prev.pvd

#save last dump into restart
dump_files=($(ls -1v dump-*))
N_dumps=${#dump_files[@]}
if (( $N_dumps > 0)); then
	last_file=${dump_files[${N_dumps} - 1]}
fi
if (( $N_dumps > 0 )); then
  echo "copying ${last_file} to restart ..."
  cp  ${last_file} restart
fi

#run the code
if [[ "$mode" == "qsub" ]] ; then
  mpirun ./a.out ${args} >>log 2>&1
else
  srun ./a.out ${args} >>log 2>&1
fi
