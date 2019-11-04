CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm $(OPENGLIBS)
my_embed:
	rm -rf .qcc* br* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) my_embed.c $(CINCL) $(CLIBS)
	./a.out > out 2> log
