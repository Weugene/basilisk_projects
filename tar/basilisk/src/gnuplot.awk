# parser for gnuplot commands in Literate C pages

BEGIN { nplots = 0 }

# This works around a bug in gnuplot 5, where font units are missing
function fixfonts(output)
{
    print "! sed -i 's/font-size=\"\\([0-9.]*\\)\"/font-size=\"\\1pt\"/g' " \
	output;
}

function defaults()
{
    printf "set pointsize 0.75; ";
}

/^[ \t]*~~~/ {
    if (gnuplot) {
	gnuplot = 0;
	print "set output";
	if (match (output, ".*\\.png"))
	    print "! mogrify -trim" output;
	else if (match (output, ".*\\.svg"))
	    fixfonts(output);
	else if (output == "")
	    fixfonts("_plot" nplots++ ".svg");
    }
}

{
    if (gnuplot) {
	if (match($0, "^[ \t]*reset")) {
	    printf "reset; ";
#	    printf "load '~/.gnuplot'; ";
	    defaults();
	}
	else {
	    if (match($0, "^[ \t]*set[ \t]+output[ \t]+(.*)", a))
		output = a[1];
	    print $0;
	}
    }
    else if (!match($0, "^[ \t]*~~~"))
	print "#", $0;
}

/^[ \t]*~~~gnuplot/ {
    gnuplot = 1; output = "";
    printf "set output '_plot" nplots ".svg'; ";
    defaults();
}
