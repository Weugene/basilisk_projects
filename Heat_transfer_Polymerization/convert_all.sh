#!/bin/bash
set -x
set -o nounset
echo "Usage:"
#parameters for solver
script=./convert_heat.sh
exec=./convertdump2pvd
list=( $(ls -d Tcyl=*) );
length=${#list[@]}

for dir in "${list[@]}"; do
  echo $dir;
  (cd $dir && cp ../$exec . &&
   ln -sf ../${script} . &&
   ./${script}
   )&
done