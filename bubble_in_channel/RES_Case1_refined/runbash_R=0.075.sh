#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
level=11
Re=0.5 #it will be overwritten
#Rad=0.3
#plist=(0.93 1.85 3.7 7.4)

#Rad=0.15
#plist=(0.74 1.85 3.7 9.3)

Rad=0.075
plist=(0.3 1.48 2.96 14.8)
for z in "${plist[@]}"; do
	echo "./a.out arg=${z} ${Rad} ${level}"
	dir="res_${Rad}_${z}_{$level}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${Rad}" "${z}" "${Re}" "${level}" > "out_${Rad}_${z}_${level}" 2> "log_${Rad}_${z}_${level}" &
	)
done
