#!/bin/bash
if [ $# -eq 0 ]; then
    echo "Usage: $0 level"
else
    level=$1
    solver=./a.out
    plist=(1 1.5 2 3 4)
    dlist=(0.0005 0.0008 0.001 0.0012)

    for b in "${dlist[@]}"; do
        for a in "${plist[@]}"; do
            dir="res_${a}_${b}_${level}"
            echo "dir=$dir: ./a.out args=$a $b"
            mkdir "$dir" || continue
            cp "$solver" "$dir"
            (
                cd "$dir" || exit
                $solver "$a" "$b" "$level" > "out_${a}_${b}_${level}" 2> "log_${a}_${b}_${level}" &
            )
        done
    done
fi
