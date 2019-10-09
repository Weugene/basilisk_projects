CFILE=three_phase_test_curved_solid
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)

three_phase_test_curved_solid:
	rm -rf enc* log .qcc*
	qcc  $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out
