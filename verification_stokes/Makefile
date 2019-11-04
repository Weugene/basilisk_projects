CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm $(OPENGLIBS)
verification_stokes:
	rm -rf .qcc* br* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) verification_stokes.c $(CINCL) $(CLIBS)
	./a.out > out 2> log
