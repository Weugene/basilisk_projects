#!/bin/bash
echo "Usage: $0 Rmu. Total input args $#"
solver=./a.out
level=10
Rmu=0.01
Re=0.1 #it will be overwritten
Rad=0.45
plist=(0.005 0.01 0.015 0.02 0.035 0.045 0.1 0.25)
for ca in "${plist[@]}"; do
	echo "./a.out arg=${Rad} ${ca} ${Re} ${level}" "${Rmu}"
	dir="res_${Rad}_${ca}_${Re}_${level}_${Rmu}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${Rad}" "${ca}" "${Re}" "${level}" "${Rmu}"> "out_${Rad}_${ca}_${Re}_${level}_${Rmu}" 2> "log_${Rad}_${ca}_${Re}_${level}_${Rmu}" &
	)
done
