#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
plist=(7 8 9 10 11)

for a in "${plist[@]}"; do
	echo "./a.out arg=$a"
	dir="res_l_$a"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "$a"  > "out_$a" 2> "log_$a" &
	)
done
