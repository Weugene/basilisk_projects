#!/bin/sh

showfiles()
{
    if ! darcs show files > /dev/null 2>&1; then
	ls *.$1 2> /dev/null | sed 's/\(.*\)/	\1 \\/g'
    else
	ROOT=`darcs show repo | grep Root | awk '{print $2}'`
	DIR=`echo $PWD | sed "s|$ROOT|.|"`
	darcs show files | grep $DIR'/[^/]*\.'$1'$' | \
	    sed 's|'$DIR'/\(.*\)|	\1 \\|g'
    fi
}

# join lines delimited by \\n characters
singleline()
{
    sed -e ':a
N
$!ba
s/\\\n//g
'
}

echo "updating Makefile.tests"
# create symbolic links between pages and sources if necessary
for f in `showfiles '[ch]\.page'` \\; do
    if test "$f" != "\\"; then
	short=`echo $f | sed 's/\.\([ch]\)\.page/\.\1/'`
	if test ! -f $short; then
	    ln -s $f $short
	fi
    fi
done
(
    echo "# Automatically generated using 'make Makefile.tests'"
    echo "# DO NOT EDIT, edit 'Makefile' instead"
    DIR=`basename $PWD`
    echo "ALLTESTS = \\"
    (showfiles c
     showfiles 'c\.page' | sed 's/\.page / /g'
     grep '^[a-zA-Z_0-9-]*\.*tst[ ]*:' Makefile | \
	 sed -n 's/\(^[a-zA-Z_0-9-]*\)\.*tst[ ]*:.*/	\1.c \\/p'
    ) | sort | uniq
    echo ""
    echo "plots: \\"
    showfiles plot | sed 's/\(.*\)\.plot/\1\/plot.png/g'
    showfiles 'c.page' | sed 's/\(.*\)\.c\.page/\1\/plots/g'
    echo ""
    echo "TESTS = \\"
    singleline < Makefile | 		                   \
	grep '^check:' | tr ' ' '\n' |                     \
	sed -n 's/[ 	]*\([a-zA-Z_0-9-]*\..tst\)[ 	]*/\1/p' |  \
	sed 's/\(.*\)\..tst/	\1.c \\/g'
    singleline < Makefile | 		                   \
	grep -v '^check:' | grep '^[^.]*:.*' | tr ' ' '\n' |     \
	sed -n 's/[ 	]*\([a-zA-Z_0-9-]*\.tst\)[ 	]*/\1/p' |    \
	sed 's/\(.*\)\.tst/	\1.c \\/g'
    echo ""
    echo "SPECIAL_TESTS = \\"
    sed -n 's/.*:[ 	]*\([a-zA-Z_0-9-]*\.ctst\)/\1/p' Makefile | \
	sort | uniq | sed 's/\(.*\)/	\1 \\/g'
    echo ""
    echo "ALLPAGES = \\"
    showfiles 'page'
    echo ""
) > Makefile.tests

rm -f Makefile.deps
