#!/bin/bash
#set -x
ulimit -n 10000
options="-y -framerate 10 -c:v libx264 -pix_fmt yuv420p -an -fs 10M "

[[ $# < 1 ]] && { echo "Usage: ./$(basename $0) <input*.format>  [<crf>]"; exit 1; }

echo ${options}
infn="$1"
extension="${infn##*.}"
fn="${infn%.*}"
pattern="$(echo "$fn"| sed 's/\*//')"
echo "infn=$infn ext=$extension filename=$fn pattern=$pattern"
#acceleration
va=2.0

exctractnum(){
	echo "$1" | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'| grep -m1 ""
}

list_files=($(ls $infn | sort -n -t = -k 2 ))
echo "list files:" ${list_files[@]}
mkdir -p fig

list_timetmp=()
for i in "${!list_files[@]}"; do
	tm=$(exctractnum ${list_files[$i]})
	echo $tm
	list_timetmp+=($tm)
done

list_time=( $( printf "%s\n" "${list_timetmp[@]}" | sort -n ) )

echo ${list_timetmp[@]} ${list_time[@]}
dt_list=()
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
        dt=$(echo "(${list_time[$i]}-(${list_time[$i-1]}))/${va}" | bc -l | sed -e 's/^\./0./' -e 's/^-\./-0./')
        dt_list+=("duration $dt")
        echo $dt, ${dt_list[$i]}
    fi
done

for i in "${!list_files[@]}"; do
	nm=${list_files[$i]}
	time=$(exctractnum $nm)
	filename="${nm%.*}"
	echo "***" $nm $filename.png
	#convert -density 600 $nm fig/$filename.png
	convert -density 600 $nm  fig/$filename.png
	#convert -density 600 $nm -gravity North  -pointsize 10 -annotate +0+100 "t=$time" fig/$filename.png   #-font  Times-Roman
done


rm -rf time_files.txt
list_file=($(ls fig/${fn}*.png | sort -n -t = -k 2 ))
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
		echo "${dt_list[$i-1]}"
        echo ${dt_list[$i-1]} >> time_files.txt
    fi
    echo "file '${list_file[$i]}'" >> time_files.txt
done

ffmpeg  -f concat -safe 0 -i time_files.txt ${options} -vf "scale=1920:1440:force_original_aspect_ratio=decrease:eval=frame,pad=1920:1440:-1:-1:color=white" ${pattern}.mp4;
