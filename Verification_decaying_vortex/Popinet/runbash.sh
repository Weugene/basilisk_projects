#!/bin/bash
echo "Usage: No arguments"
solver=./popinet.out
#levels=(4)
levels=(4 5 6 7 )
echo "OPENGLIBS=$OPENGLIBS"
for a in "${levels[@]}"; do
	echo "./a.out arg=$a"
	dir="res_l_${a}"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "${a}"  >  "out_${a}" 2> "log_${a}" &
	)
done
