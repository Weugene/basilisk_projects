#set -x
echo "$0 prefix root"
ulimit -n 10000
options="-y -framerate 10 -c:v libx264 -pix_fmt yuv420p -an -fs 10M "

echo ${options}
prefix=$1
root=$2
scalex=$3
scaley=$4
scale_opt=" -vf scale=${scalex}:${scaley}:force_original_aspect_ratio=decrease:eval=frame,pad=${scalex}:${scaley}:-1:-1:color=white"
if [[ -z $scalex ||  -z $scaley ]]; then echo "scalex is not defined, then set default"; scale_opt=''; fi

dn="${prefix}_t=*${root}.png"
#pic_t=4.59881ux.png
echo "SEE:" $(ls $dn)
va=2.0
list_time=($(ls $dn | sort -n -t = -k 2 | sed "s/${prefix}_t=//" | sed "s/${root}.png//"));
echo ${list_time[@]}
dt_list=()
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
        dt=$(echo "(${list_time[$i]}-(${list_time[$i-1]}))/${va}" | bc -l | sed -e 's/^\./0./' -e 's/^-\./-0./')
        dt_list+=("duration $dt")
        echo $dt, $dt_list[$i]
    fi
done

mkdir -p fig

echo $root;
rm ${root}_files.txt
list_file=($(ls ${dn} | sort -n -t = -k 2 ))
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
        echo ${dt_list[$i-1]} >> ${root}_files.txt
    fi
    echo "file '${list_file[$i]}'" >> ${root}_files.txt
    #nm=${list_files[$i]}
    #time=${list_time[$i]}
    #filename="${nm%.*}"
    #convert -density 600 $nm -gravity North  -pointsize 10 -annotate +0+100 "t=$time" fig/$filename.png
done


ffmpeg  -f concat -safe 0 -i ${root}_files.txt ${options} ${scale_opt} video_${root}.mp4;

