#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
levels=(6 7 8 9 )
eta_s=(1e-3 1e-4 1e-5 1e-6)
dtlimiter=(0 1)
echo "OPENGLIBS=$OPENGLIBS"
for dtlim in "${dtlimiter[@]}"; do
for eta in "${eta_s[@]}"; do
for a in "${levels[@]}"; do
	echo "./a.out arg=$a ${eta} ${dtlim}"
	dir="res_l_${a}_${eta}_${dtlim}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${a}" "${eta}" "${dtlim}" > "out_${a}_${eta}_${dtlim}" 2> "log_${a}_${eta}_${dtlim}" &
	)
done
done
done
