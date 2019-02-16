# Processes $foo()$ include commands in .st files

BEGIN {
    path = "/usr/share/gitit/data/templates/"
}

{
    if ($0 ~ /\$messages:listitem\(\)\$/) {
	# This is a strange command we don't understand, we just remove it
    }
    else if ($0 ~ /^[ \t]*\$[a-zA-Z0-9_]+\(\)\$[ \t]*$/) {
	gsub ("^[ \t]*\\$", "")
	gsub ("\\(\\)\\$[ \t]*$", "")
	static = $0 ".static"
	file = $0 ".st"
	system ("(test -f " static " && cat " static ") || "		\
	        "(test -f " file " && cat " file ") || "		\
		"(test -f " path file " && cat " path file ")");
    }
    else if ($0 ~ /^[ \t]*\$content\$[ \t]*$/)
	system ("cat body.static");
    else if ($0 ~ /^[ \t]*\$tabs\$[ \t]*$/)
	system ("cat tabs.static");
    else
	print $0;
}
