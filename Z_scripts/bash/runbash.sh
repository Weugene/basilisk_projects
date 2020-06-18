#!/bin/bash
echo "Usage: No arguments"
solver=./a.out
plist=(5 6 7 8 9)

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
