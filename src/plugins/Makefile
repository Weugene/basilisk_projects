all: LiterateC.o CodeBlock.o

LiterateC.o: LiterateC.hs literate-c.o
	ghc -Wall LiterateC.hs
	mv -f LiterateC.o LiterateC1.o
	ld -r literate-c.o LiterateC1.o -o LiterateC.o
	rm -f LiterateC1.o

CodeBlock.o: CodeBlock.hs codeblock.o
	ghc -Wall CodeBlock.hs
	mv -f CodeBlock.o CodeBlock1.o
	ld -r codeblock.o CodeBlock1.o -o CodeBlock.o
	rm -f CodeBlock1.o

literate-c: literate-c.c
	$(CC) -g -O2 -DSTANDALONE=1 literate-c.c -o literate-c

codeblock: codeblock.c
	$(CC) -g -O2 -DSTANDALONE=1 codeblock.c -o codeblock

codeblock-test: codeblock.c
	$(CC) -g -O2 -DTEST=1 codeblock.c -o codeblock-test

clean:
	rm -f *.o

%.o: %.c
	$(CC) -g -O2 -Wall -Wno-unused -c $<

%.o: %.hs
	ghc -g -Wall $<

%.c: %.lex
	flex $<
	@sed -e 's/^[ \t]*int \(yy[a-z_]*[ \t]*(\)/static int \1/g' \
	    -e 's/^[ \t]*YY_BUFFER_STATE \(yy[a-z_]*[ \t]*(\)/static YY_BUFFER_STATE \1/g' \
	    -e 's/^[ \t]*YY_EXTRA_TYPE \(yy[a-z_]*[ \t]*(\)/static YY_EXTRA_TYPE \1/g' \
	    -e 's/^[ \t]*FILE \*\(yy[a-z_]*[ \t]*(\)/static FILE \*\1/g' \
	    -e 's/^[ \t]*char \*\(yy[a-z_]*[ \t]*(\)/static char \*\1/g' \
	    -e 's/^[ \t]*void \*\(yy[a-z_]*[ \t]*(\)/static void \*\1/g' \
	    -e 's/extern int yylex/static int yylex/g' \
            -e 's/^[ \t]*void \(yy[a-z_]*[ \t]*(\)/static void \1/g' < lex.yy.c > $@
