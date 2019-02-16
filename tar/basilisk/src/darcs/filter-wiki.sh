#!/bin/bash
# unpull patches changing files with a .page extension or file in directories
# other than ./src/

unpull() {
    last=`echo "$DARCS_PATCHES_XML" | grep -c '<patch '`
    echo "filter-wiki.sh: unpulling the following $last 'wiki patches':\n"
    echo "$DARCS_PATCHES"
    darcs unpull -a --last=$last
    exit 1;
}

for file in $DARCS_FILES; do
    case "$file" in
	*.page)  unpull;;
	./src/*)       ;;
	*)       unpull;;
    esac
done
