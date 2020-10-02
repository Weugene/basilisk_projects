#set -x
ulimit -n 10000
options="-y -framerate 10 -c:v libx264 -pix_fmt yuv420p -an -fs 10M"
echo ${options}
dn='lambda2_t=*ux.png'
sdn='ux'
list_time=($(ls $dn | sort -n -t = -k 2 | sed "s/lambda2_t=//" | sed "s/${sdn}.png//"));
echo ${list_time[@]}
dt_list=()
for i in "${!list_time[@]}"; do
    if [[ $i > 0 ]]; then
        dt=$(echo "${list_time[$i]}-(${list_time[$i-1]})" | bc -l | sed -e 's/^\./0./' -e 's/^-\./-0./')
        dt_list+=("duration $dt")
        echo $dt, $dt_list[$i]
    fi
done
for f in Lambda2_in_bubble ux Lambda2 tracer uxSlice uxSide; do
    echo $f;
    rm ${f}_files.txt
    list_file=($(ls *$f.png | sort -n -t = -k 2 ))
    for i in "${!list_time[@]}"; do
        if [[ $i > 0 ]]; then
            echo ${dt_list[$i]} >> ${f}_files.txt
        fi
        echo "file '${list_file[$i]}'" >> ${f}_files.txt
    done

    ffmpeg  -f concat -safe 0 -i ${f}_files.txt ${options} ${f}.mp4;
    #ffmpeg -i "concat:$list" ${options} ${f}.mp4;
done;

for f in uxSlice omegaSlice; do
    echo $f;
    rm ${f}_files.txt
    list_file=($(echo $(ls lambda2_t=*${f}.png| sort -n -t = -k 2 -k 4) | sed -e "s/ /|/g"))
    for i in "${!list_time[@]}"; do
        if [[ $i > 0 ]]; then
            echo ${dt_list[$i]} >> ${f}_files.txt
        fi
        echo "file '${list_file[$i]}'" >> ${f}_files.txt
    done
    ffmpeg  -f concat -safe 0 -i ${f}_files.txt ${options} ${f}.mp4;
done;

#echo "M_bubble_t=...";
#list=$(echo $(ls M_bubble_t=*.pdf| sort -n -t = -k 2) | sed -e "s/ /|/g")
#echo "M_bubble: $list"
#ffmpeg -framerate 5 -i "concat:$list" -c:v ffv1 M_bubble.avi;

#videoConverter.sh -y -q=1 -fs 10M  *.avi;
#rm *.avi;

