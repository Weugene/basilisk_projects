ulimit -n 10000
options="-y -framerate 10 -i "concat:$list" -c:v libx264 -pix_fmt yuv420p -an -fs 10M"
for f in {"ux","Lambda2","tracer","uxSlice","uxSide"}; do
    echo $f;
    list=$(echo $(ls lambda2_t=*${f}.png| sort -n -t = -k 2) | sed -e "s/ /|/g")
    echo "$f: $list"
    ffmpeg ${options} ${f}.mp4;
done;

for f in {"uxSlice","omegaSlice"}; do
    echo $f;
    list=$(echo $(ls lambda2_t=*${f}.png| sort -n -t = -k 2 -k 3) | sed -e "s/ /|/g")
    echo "$f: $list"
    ffmpeg ${options} ${f}.mp4;
done;

#echo "M_bubble_t=...";
#list=$(echo $(ls M_bubble_t=*.pdf| sort -n -t = -k 2) | sed -e "s/ /|/g")
#echo "M_bubble: $list"
#ffmpeg -framerate 5 -i "concat:$list" -c:v ffv1 M_bubble.avi;

#videoConverter.sh -y -q=1 -fs 10M  *.avi;
#rm *.avi;

