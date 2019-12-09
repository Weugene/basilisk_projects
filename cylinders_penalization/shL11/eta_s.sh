#!/bin/bash -v

solver=./a.out
level=10
plist=(1e-6 1e-5 1e-4 1e-3 1e-2)

for eta in "${plist[@]}"; do
	echo "eta=$eta level=$level"
	dir="res$eta"
	mkdir "$dir" || continue
	cp "$solver" "$dir"
	(
		cd "$dir" || exit
		$solver "$eta" $level > "out_$eta" 2> "log_$eta" &
	)
done
wait
for eta in "${plist[@]}"; do
	rm "./res$eta/$solver"
done

