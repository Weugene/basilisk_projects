#!/bin/bash
echo "from to step first_name";
for n in $(seq -f "%04g" $1 $3 $2); do rsync -av l:~/wbasilisk/${PWD##*/}/$4*${n}* .; done
