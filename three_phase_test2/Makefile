CFILE=three_phase_test2
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)

three_phase_test2:
	rm -rf enc* log .qcc* dump
	qcc  $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out
