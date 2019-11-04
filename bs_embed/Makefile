CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm $(OPENGLIBS)
basilisk_embed:
	rm -rf .qcc* br* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) basilisk_embed.c $(CINCL) $(CLIBS)
	./a.out
