#!/bin/bash
echo "Usage: $0 Rrho Rmu maxlevel. Total input args $#"
source $HOME/.bashrc
solver=./a.out
Rrho=$1
Rmu=$2
maxlevel=$3
Radlist=(1)
Relist=(0.1 1 10 100)
Calist=(0.01 0.1 1)
for rad in "${Radlist[@]}";do
for re in "${Relist[@]}";do
for ca in "${Calist[@]}"; do
	echo "./a.out arg=${re} ${ca} ${rad} ${Rrho}" "${Rmu}" "${maxlevel}"
	var=${re}_ca=${ca}_r=${rad}_Rrho=${Rrho}_Rmu=${Rmu}_l=${maxlevel}
   	dir="res_$var"
   	echo "${dir}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${re}" "${ca}" "${rad}" "${Rrho}" "${Rmu}" "${maxlevel}"> "out_${var}" 2> "log_${var}" &
	)
done
done
done
wait
