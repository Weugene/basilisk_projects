#set -x
ulimit -n 10000
options="-y -framerate 10 -c:v libx264 -pix_fmt yuv420p -an -fs 10M "

echo ${options}
sdn='_uxSide'
dn="lambda2_t=*${sdn}.png"

va=2.0
list_time=($(ls $dn | sort -n -t = -k 2 | sed "s/lambda2_t=//" | sed "s/${sdn}.png//"));
echo ${list_time[@]}
dt_list=()
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
        dt=$(echo "(${list_time[$i]}-(${list_time[$i-1]}))/${va}" | bc -l | sed -e 's/^\./0./' -e 's/^-\./-0./')
        dt_list+=("duration $dt")
        echo $dt, $dt_list[$i]
    fi
done
for f in uxSide ; do
    echo $f;
    rm ${f}_files.txt
    list_file=($(ls *$f.png | sort -n -t = -k 2 ))
    for i in "${!list_time[@]}"; do
        if [[ $i > 0 ]]; then
            echo ${dt_list[$i-1]} >> ${f}_files.txt
        fi
        echo "file '${list_file[$i]}'" >> ${f}_files.txt
    done

    ffmpeg  -f concat -safe 0 -i ${f}_files.txt ${options} -vf "scale=3108:1168:force_original_aspect_ratio=decrease:eval=frame,pad=3180:1168:-1:-1:color=white" ${f}.mp4;
done;

for f in uxSlice omegaSlice; do
    echo $f;
    rm ${f}_files.txt
    list_file=($(ls lambda2_t=*${f}.png| sort -n -t = -k 2 -k 4 | sed -e "s/ /|/g"))
    for i in "${!list_time[@]}"; do
        if [[ $i > 0 ]]; then
            echo ${dt_list[$i-1]} >> ${f}_files.txt
        fi
        echo "file '${list_file[$i]}'" >> ${f}_files.txt
    done;
    ffmpeg  -f concat -safe 0 -i ${f}_files.txt ${options} ${f}.mp4;
done;
