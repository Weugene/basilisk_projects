#!/bin/bash
set -x
set -o nounset
list=( $(ls dump-* | sort -n -t - -k 2) );
length=${#list[@]}


pvd="dump2pvd.pvd"
#subn="${pvd%.*}"

subn_pvd="dump2pvd_heat"
pvdconvertD2P="$subn_pvd.pvd"
time_of_debug=true

if $time_of_debug; then
    rm $pvd $pvdconvertD2P;
    text1='<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
    <Collection>'
    text3='    </Collection>
</VTKFile>'

    for (( i = 0; i < length; i+=1 )); do
      echo "i=$i";
      dump=${list[i]}
      t=$(echo "$dump" | sed 's/dump-//')
      echo "$dump at time=$t";
      #sleep 4
      id_out=$i
      ./convertdump2pvd $dump $i 14
      text2_pvd=""
      for (( j = 0; j <= i; j+=1 )); do
        dumpj=${list[j]}
        tj=$(echo "$dumpj" | sed 's/dump-//')
        text2_pvd="$text2_pvd   <DataSet timestep=\"${tj}\" part=\"0\" file=\"res/${subn_pvd}_0_$(printf %04d ${j}).pvtu\"/>"
      done;
      echo $text1 >> $pvd
      echo $text2_pvd >> $pvd
      echo $text3 >> $pvd

      # Wait for all jobs to complete
      jobs
      wait
      jobs
#      sleep 5
      mv $pvd $pvdconvertD2P
    done

fi
