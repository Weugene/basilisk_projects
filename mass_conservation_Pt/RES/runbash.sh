#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
plist=(0.08 0.125 0.2)

for a in "${plist[@]}"; do
	echo "./a.out arg=$a"
	dir="res_$a"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "$a"  > "out_$a" 2> "log_$a" &
	)
done
