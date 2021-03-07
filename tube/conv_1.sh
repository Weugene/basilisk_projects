#!/bin/bash
set -x
set -o nounset
echo "NOTE: reads only dump2pvd_compressed file!!!"
echo "Usage ./conv.sh Ncase maxlevel volumetric_repr"
bubcase=$1
maxlevel=$2
volumetric_repr=$3
echo "in: Ncase=$bubcase maxlevel=$maxlevel"
for f in convert_single extract_iso_volume.py; do
    test -e ${f} && echo "${f} exists" || echo "${f} does not exist. exit" || exit;
done

list=( $(ls dump-* | sort -n -t - -k 2) );
length=${#list[@]}

DISPLAY_DEF=17
d=$(ls -lha /tmp/.X11-unix/ |grep ${USER} | awk '{print $9}' | sed 's/X//' | sed 's/=//'| tail -n 1)
[[ -z "$d" || "$USER" == 'weugene' ]] && echo "DISPLAY not set!" && vncserver :${DISPLAY_DEF} &&  d=${DISPLAY_DEF}
echo "DISPLAY=:${d}"

pvd="dump2pvd.pvd"
#subn="${pvd%.*}"
subn_isovolume="save_isovolume"
subn_pvd="dump2pvd_compressed"
isopvd="$subn_isovolume.pvd"
pvdconvertD2P="$subn_pvd.pvd"
time_of_debug=true

if $time_of_debug; then
    rm $pvd;
    text1='<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
    <Collection>'
    text3='    </Collection>
</VTKFile>'
    echo $text1 >> $isopvd
    for (( i = 0; i < length; i+=1 )); do
      echo "i=$i";
      dump=${list[i]}
      t=$(echo "$dump" | sed 's/dump-//')
      echo "$dump at time=$t";
      #sleep 4
      ./convert_single $dump $i $bubcase $maxlevel
      text2_pvd="        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn_pvd}_0_$(printf %04d $i).pvtu\"/>" >> $pvd      
      text2_isovolume="        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn_isovolume}_0_$(printf %04d $i).pvtu\"/>" >> $pvd
      echo $text1 >> $pvd
      echo $text2_pvd >> $pvd
      echo $text3 >> $pvd
      echo $text2_isovolume >> $isopvd
      # Wait for all jobs to complete
      jobs
      wait
      jobs
      sleep 5
      mv $pvd $pvdconvertD2P
      DISPLAY=:${d} VTK_USE_LEGACY_DEPTH_PEELING=1 pvpython extract_iso_volume.py -infn $pvdconvertD2P -outfn $subn_isovolume -maxlevel $maxlevel -iter $i -volumetric_repr $volumetric_repr
    done
    echo $text3 >> $isopvd
fi
