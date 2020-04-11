#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
level=10

Re=0.1
Rad=0.3
plist=(0.0005 0.001 0.005 0.01 0.05 0.1)
for a in "${plist[@]}"; do
	echo "./a.out arg=${Rad} ${a} ${Re} ${level}"
	dir="res_${Rad}_${a}_${Re}_${level}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${Rad}" "${a}" "${Re}" "${level}" > "out_${Rad}_${a}_${Re}_${level}" 2> "log_${Rad}_${a}_${Re}_${level}" &
	)
done
