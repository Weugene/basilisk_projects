pattern=$1
ext="${pattern##*.}"
fn="${pattern%.*}"
for n in $(ls $pattern); do 
	name=$(basename $n .$ext); 
	echo $n $name; 
	convert -density 500  ${name}.${ext}  -resize 1024x1024 ${name}.png; 
done
