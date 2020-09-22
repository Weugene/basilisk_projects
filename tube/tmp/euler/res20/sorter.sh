
for ff in { uxSlice omegaSlice }; do
    for i in $(ls *${ff}.png); do
        echo "$i" | sed 's/:*lambda2_t=//' | sed 's/:*_n=/ /' | sed "s/:*_${ff}.png//"; >> "${ff}".txt;
    done
done
