all: literate-c codeblock pagemagic sanitize urldecode

CFLAGS += -Wall -Wno-unused -g -O2

codeblock: codeblock.c

literate-c: literate-c.c

# add flex -d below for debugging
%.c: %.lex Makefile
	flex $<
	@sed -e 's/^[ \t]*int \(yy[a-z_]*[ \t]*(\)/static int \1/g' \
	    -e 's/^[ \t]*YY_BUFFER_STATE \(yy[a-z_]*[ \t]*(\)/static YY_BUFFER_STATE \1/g' \
	    -e 's/^[ \t]*YY_EXTRA_TYPE \(yy[a-z_]*[ \t]*(\)/static YY_EXTRA_TYPE \1/g' \
	    -e 's/^[ \t]*FILE \*\(yy[a-z_]*[ \t]*(\)/static FILE \*\1/g' \
	    -e 's/^[ \t]*char \*\(yy[a-z_]*[ \t]*(\)/static char \*\1/g' \
	    -e 's/^[ \t]*void \*\(yy[a-z_]*[ \t]*(\)/static void \*\1/g' \
	    -e 's/extern int yylex/static int yylex/g' \
            -e 's/^[ \t]*void \(yy[a-z_]*[ \t]*(\)/static void \1/g' < lex.yy.c > $@
