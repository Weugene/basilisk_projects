set -o nounset
echo "Usage $0 outfilename deltafile. Be careful! Avoid usage copy in deltafile"
echo $1, $2
file=$1
deltafile=$2

#time1=($(grep "timestep" $file | awk '{print $2}'))
time1=($(grep -e 'slice n=0' -e 'times' $file  | grep -e 'slice' -B 1 | sed '/^--$/d' | grep timesteps | sed 's/timesteps= \[//' |sed 's/\]//'))
volume=($(grep -v "cand" $file | grep  -e "refined_x_mean" | awk '{print $8}'))
umean=($(grep -v "cand" $file | grep  -e "refined_x_mean" | awk '{print $4}' | sed 's/\[//'| sed 's/,//'))
xmean=($(grep -v "cand" $file | grep  -e "refined_x_mean" | awk '{print $2}'))
xpeak=($(grep -v "cand" $file | grep -e "x_peak" | awk '{print $2}'))
ypeak=($(grep -v "cand" $file | grep -e "x_peak" | awk '{print $4}'))

time2=($(cat $deltafile | awk '{print $1}'))
delta_min=($(cat $deltafile | awk '{print $2}'))
delta_mean=($(cat $deltafile | awk '{print $3}'))
delta_max=($(cat $deltafile | awk '{print $4}'))
delta_min_smooth=($(cat $deltafile | awk '{print $5}'))
delta_max_smooth=($(cat $deltafile | awk '{print $6}'))
x0=($(cat $deltafile | awk '{print $7}'))
y0=($(cat $deltafile | awk '{print $8}'))
xN=($(cat $deltafile | awk '{print $9}'))
yN=($(cat $deltafile | awk '{print $10}'))

echo "t	x_tail	x_peak	y_peak	x_mean	x_nose	x_nose_ISC	volume	UmeanV	delta_min	delta_mean	delta_max	delta_min_smooth	delta_max_smooth"
for j in "${!time2[@]}"; do
    for i in "${!time1[@]}"; do
		if  (( $(echo "${time1[$i]} == ${time2[$j]}" |bc -l) )); then
		   echo "${time1[$i]}	${x0[$j]}	${xpeak[$i]}	${ypeak[$i]}	${xmean[$i]}	${xN[$j]}	?	${volume[$i]}	${umean[$i]}	${delta_min[$j]}	${delta_mean[$j]}	${delta_max[$j]}	${delta_min_smooth[$j]}	${delta_max_smooth[$j]}";
		fi
	done
done

#t	x_tail	x_peak	y_peak	x_mean	x_nose	x_nose_ISC	volume	UmeanV	delta_min	delta_mean	delta_max	delta_min_smooth	delta_max_smooth
