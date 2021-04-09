#!/bin/bash
set -x
set -o nounset
echo "Usage: Tcyl Tin maxlevel dr dx dy"

Tcyl=$1
Tin=$2
maxlevel=$3
iter_fp=0
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
mkdir ${tmp} || exit
cd ${tmp} || exit
cp ../a.out .
args="${Tcyl} ${Tin} ${maxlevel} ${iter_fp} ${ratio_Rbmin} ${ratio_Rbmax} ${ratio_dist_x} ${ratio_dist_y} ${ratio_front_x} \
${cyl_x} ${Nb} ${Ncx} ${Ncy} ${TOLERANCE_P} ${TOLERANCE_V} ${TOLERANCE_T} \
${Htr} ${Arrhenius_const} ${Ea_by_R} ${dev_r} ${dev_x} ${dev_y}"

echo "args:${args}"

echo "BASILISK=$BASILISK"

if [[ "$mode" == "qsub" ]] ; then
  mpirun ./a.out ${args} > out 2> log
else
  srun ./a.out ${args} > out 2> log
fi
