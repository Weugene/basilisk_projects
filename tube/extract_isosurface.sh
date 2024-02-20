#!/bin/bash
set -x
set -o nounset
echo "NOTE: reads only dump2pvd_compressed file!!!"
echo "Usage ./extract_isosurface.sh Ncase maxlevel iter_cur"
module load visualization/vnc-server || continue
bubcase=$1
maxlevel=$2
iter_cur=$3
range_colorbar=$4 # 0 means automatic range for colorbar
outPrefix=""
dumpPattern="dump-*"
echo "in: Ncase=$bubcase maxlevel=$maxlevel iter_cur=$iter_cur range_colorbar=$range_colorbar"
for f in convert_single extract_isosurface.py; do
    test -e ${f} && echo "${f} exists" || echo "${f} does not exist. exit" || exit;
done

list=( $(ls dump-* | sort -n -t - -k 2) );
length=${#list[@]}

DISPLAY_DEF=19
d=$(ls -lha /tmp/.X11-unix/ |grep ${USER} | awk '{print $9}' | sed 's/X//' | sed 's/=//'| tail -n 1)
[[ -z "$d" || "$USER" == 'weugene' ]] && echo "DISPLAY not set!" && vncserver :${DISPLAY_DEF} &&  d=${DISPLAY_DEF} ||  d=${DISPLAY_DEF}
echo "DISPLAY=:${d}"

pvd="dump2pvd.pvd"
#subn="${pvd%.*}"


subn_pvd="dump2pvd_compressed"
pvdconvertD2P="$subn_pvd.pvd"
save_all_data=false

for (( i = $iter_cur; i < $length; i+=1 )); do
  echo "i=$i";
  dump=${list[i]}
  t=$(echo "$dump" | sed 's/dump-//')
  echo "$dump at time=$t";
  #sleep 4
  id_out=$i
  if $save_all_data; then
      ./convert_single $dump $i $bubcase $maxlevel
  else
      ./convert_single $dump 0 $bubcase $maxlevel
      id_out=0
  fi
  # Wait for all jobs to complete
  wait
  sleep 5
  DISPLAY=:${d} VTK_USE_LEGACY_DEPTH_PEELING=1 pvpython extract_isosurface.py --infn "$pvdconvertD2P" --dumpPattern "$dumpPattern" --maxlevel $maxlevel --iter $i --rangeColorbar $range_colorbar

done