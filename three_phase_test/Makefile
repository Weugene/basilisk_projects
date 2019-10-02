CFILE=three_phase_test
CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)

$(CFILE):
	rm -rf enc* log .qcc*
	qcc  $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out