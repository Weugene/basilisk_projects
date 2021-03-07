BEGIN{ FS=", " }

/@hal{.*,.*}/{
    gsub ("}", "", $2);
    gsub (".*{", "", $1);
    system ("wget -q -O - https://hal.archives-ouvertes.fr/" $2 "/bibtex | " \
	    "sed 's/^\\(@.*\\){.*,/\\1{" $1 ",/g'");
    next;
}

{ print $0 }
