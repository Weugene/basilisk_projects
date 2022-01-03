#!/bin/bash
set -x
set -o nounset
dirs=($(ls -d ../Tcyl*))
for dir in "${dirs[@]}"; do
  #extract all numbers from the line
  arr=($(echo "${dir}" | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'))
  base_dir=$(basename $dir)
  Tcyl=${arr[0]}
  Tin=${arr[1]}
  maxl=${arr[2]}
  rdx=${arr[3]}
  rdy=${arr[4]}
  rfx=${arr[5]}
  Ncx=${arr[6]}
  Ncy=${arr[7]}
  dr=${arr[8]}
  dx=${arr[9]}
  dy=${arr[10]}
  new_fp=permiability_${base_dir}.txt
  cp ${dir}/permiability.txt ${new_fp}
  cat ${new_fp} | sed 's/vol_max/vol_max\n/' > tmp_${new_fp}
  mv tmp_${new_fp} ${new_fp}
done