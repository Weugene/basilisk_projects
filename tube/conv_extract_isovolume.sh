#!/bin/bash
set -x
set -o nounset
echo "NOTE: reads only dump2pvd_compressed file!!!"
echo "Usage ./conv_extract_isovolume.sh Ncase maxlevel volumetric_repr iter_cur"
module load visualization/vnc-server
bubcase=$1
maxlevel=$2
volumetric_repr=$3
iter_cur=$4
echo "in: Ncase=$bubcase maxlevel=$maxlevel volumetric_repr=$volumetric_repr"
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
if (($volumetric_repr==1)); then
	subn_isovolume="save_isovolume"
	format_out="pvtu"
else
	subn_isovolume="save_isosurface"
	format_out="vtp"
fi
subn_pvd="dump2pvd_compressed"
isopvd="$subn_isovolume.pvd"
isopvdtail="${subn_isovolume}_tail.pvd"
pvdconvertD2P="$subn_pvd.pvd"
time_of_debug=true
save_all_data=false
if $time_of_debug; then
    rm $pvd $isopvd $isopvdtail $pvdconvertD2P;
    text1='<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
    <Collection>'
    text3='    </Collection>
</VTKFile>'
    
    for (( i = $iter_cur; i < length; i+=1 )); do
      echo "i=$i";
      dump=${list[i]}
      t=$(echo "$dump" | sed 's/dump-//')
      echo "$dump at time=$t";
      #sleep 4
      id_out=$i
      if $save_all_data; then
          ./convert_single $dump $i $bubcase $maxlevel
      else
          ./convert_single $dump 0 $bubcase $maxlevel	  
	  id_out=0
      fi
      
      text2_pvd="        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn_pvd}_0_$(printf %04d ${id_out}).pvtu\"/>" >> $pvd      
      
      rm $pvd
      echo $text1 >> $pvd
      echo $text2_pvd >> $pvd
      echo $text3 >> $pvd
      
      # Wait for all jobs to complete
      jobs
      wait
      jobs
      sleep 5
      mv $pvd $pvdconvertD2P
      DISPLAY=:${d} VTK_USE_LEGACY_DEPTH_PEELING=1 pvpython extract_iso_volume.py -infn $pvdconvertD2P -outfn $subn_isovolume -maxlevel $maxlevel -iter $i -volumetric_repr $volumetric_repr
      rm $isopvd $isopvdtail
      echo $text1 >> $isopvd
      for (( j = 0; j <= i; j+=1 )); do
          dump=${list[j]}
          t=$(echo "$dump" | sed 's/dump-//')
          text2_isovolume="        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn_isovolume}_0_$(printf %04d $j).${format_out}\"/>" >> $pvd
          echo $text2_isovolume >> $isopvd
      done
	  echo $text3 >> $isopvd

	  # bubble surface tail
	  echo $text1 >> $isopvdtail
      for (( j = 0; j <= i; j+=1 )); do
          dump=${list[j]}
          t=$(echo "$dump" | sed 's/dump-//')
          text2_isovolume="        <DataSet timestep=\"${t}\" part=\"0\" file=\"res/${subn_isovolume}_tail_0_$(printf %04d $j).${format_out}\"/>" >> $pvd
          echo $text2_isovolume >> $isopvdtail
      done
      echo $text3 >> $isopvdtail
    done
    
fi
