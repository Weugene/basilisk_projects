src=$1
grep "incl .*/"$src \
    $BASILISK/*.tags \
    $BASILISK/navier-stokes/*.tags \
    $BASILISK/examples/*.tags \
    $BASILISK/test/*.tags | \
    awk -v basilisk=$BASILISK '
            function title(fname) {
              while ($1 != "#") {
                if (getline <fname == 0)
                  return fname;
              }
              gsub("# ", "", $0);
              gsub("*/$", "", $0);
              gsub("\"\"\"$", "", $0);
              return $0;
            }
            {
              gsub(basilisk, "", $1);
              gsub(".tags:.*", "", $1);
              used = "/src" $1;
              lineno = $4;
              print "used " title(basilisk $1 ".page") "\t" used " " lineno;
            }'
