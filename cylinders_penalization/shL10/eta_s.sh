#!/bin/bash
echo "Usage: $0 level; No args=> cp logs and outs. Total input args $#"
solver=./a.out
level=$1
plist=(1e-15 1e-6 1e-5 1e-4 1e-3 1e-2)

if [ $# -gt 0 ]; then
	for eta in "${plist[@]}"; do
		echo "eta=$eta level=$level"
		dir="res$eta"
		mkdir "$dir" || continue
		cp "$solver" "$dir"
		(
			cd "$dir" || exit
			$solver "$eta" "$level" > "out_$eta" 2> "log_$eta" &
		)
	done
else
	for eta in "${plist[@]}"; do
		dir="res$eta"
		cd "$dir" || exit
		cp log_$eta out_$eta ../
		cd ..
		#rm "./res$eta/$solver"
	done
fi
