
for f in {"uxSlice","ux","tracer","Lambda2"}; do
    echo $f.png;
    ffmpeg -framerate 5 -pattern_type glob -i "*${f}.png" -c:v ffv1 ${f}.avi;
done;
videoConverter.sh -q=1 *.avi;
rm *.avi;
