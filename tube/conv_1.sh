#!/bin/bash
set -x
set -o nounset
echo "NOTE: reads only dump2pvd_compressed file!!!"
echo "Usage ./conv.sh Nparallel Ncase shiftm shiftp maxlevel nslices"
echo "in: Nparallel=$1 Ncase=$2 shiftm=$3 shiftp=$4 maxlevel=$5 nslices=$6"
for f in convert_single lambda2_final.py png2avi.sh; do
    test -e ${f} && echo "${f} exists" || echo "${f} does not exist. exit" || exit;
done
parallel=$1
parallel=${parallel:=1}
list=( $(ls dump-* | sort -n -t - -k 2) );
length=${#list[@]}

d=$(ls -lha /tmp/.X11-unix/ |grep ${USER} | awk '{print $9}' | sed 's/X//' | sed 's/=//'| tail -n 1)
pvd="dump2pvd.pvd"
#subn="${pvd%.*}"
subn="dump2pvd_compressed"
time_of_debug=true

if $time_of_debug; then
    rm $pvd
    for (( i = 0; i < length; i+=$parallel )); do
      echo "i=$i";
      for (( j = 0; j < $parallel; j++ )); do
        echo '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
    <Collection>' >> $pvd
        k=$((i + j))
        echo "length=$length parallel=${parallel} k=$k"
        if (( k < length)); then
          dump=${list[k]}
          t=$(echo "$dump" | sed 's/dump-//')
          echo "$dump at time=$t";
          #sleep 4
          ./convert_single $dump $j $2 $3 $4 &
          echo "        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn}_0_$(printf %04d $j).pvtu\"/>" >> $pvd
        fi
      done
      # Wait for all jobs to complete
      jobs
      wait
      jobs
      echo '    </Collection>
</VTKFile>' >> $pvd
      mv "$pvd"  "$subn.pvd"
      sleep 5
      DISPLAY=:${d} VTK_USE_LEGACY_DEPTH_PEELING=1 pvpython lambda2_final.py -maxlevel $5 -nslices $6 -filename dump2pvd_compressed.pvd
    done
fi
./png2avi.sh
