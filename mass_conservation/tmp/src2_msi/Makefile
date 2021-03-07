CFLAGS += -O2

# these are not tests
EXCLUDE = qcc.c include.c rotate.c qplot.c bview-server.c

all: qcc qplot libkdt literatec draw_get.h
	@chmod +x ppm2mpeg ppm2mp4 ppm2ogv ppm2gif runtest \
		sequence bview bppm page2html
	@test -f xyz2kdt || ln -s kdt/xyz2kdt
	@test -f kdtquery || ln -s kdt/kdtquery

bview-servers: bview-server2D bview-server3D bview2D bview3D

libkdt:
	cd kdt && make

literatec:
	cd darcsit && make && cd cgi-bin && make

qcc: qcc.c include.o rotate.o config
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" \
		-DCC99="\"$(CC99)\"" \
		-DCPP99="\"$(CPP99)\"" \
		-DCADNACC="\"$(CADNACC)\"" \
		-DBASILISK="\"$(BASILISK)\"" \
		qcc.c include.o rotate.o -o qcc

rotate.c: rotate.lex
	flex -P rot -o rotate.c rotate.lex

rotate: rotate.c
	$(CC) $(CFLAGS) -DTEST rotate.c -o rotate

rotate.o: rotate.c
	$(CC) $(CFLAGS) -c rotate.c

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC) $(CFLAGS) -c qplot.c

qplot: qplot.o
	$(CC) $(CFLAGS) qplot.o -o qplot

draw_get.h: draw.h params.awk
	awk -f params.awk < draw.h > draw_get.h

bview-server2D.d: bview-server.c qcc
	qcc -MD -o $@ $<

bview-server2D: bview-server.c bview-server2D.d draw_get.h parse.h
	qcc $(CFLAGS) bview-server.c -o $@ -Lgl -lglutils $(OPENGLIBS) -lm

bview-server3D.d: bview-server.c qcc
	qcc -grid=octree -MD -o $@ $<

bview-server3D: bview-server.c bview-server3D.d draw_get.h parse.h
	qcc -grid=octree $(CFLAGS) bview-server.c -o $@ \
		-Lgl -lglutils $(OPENGLIBS) -lm

bview2D:
	ln -s bview bview2D

bview3D:
	ln -s bview bview3D

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex Makefile
	flex -o qcc.c qcc.lex

alltags:
	cd navier-stokes && make tags
	cd layered && make tags
	cd ehd && make tags
	cd examples && make tags
	cd test && make tags
	make tags
	cd navier-stokes && make itags
	cd layered && make itags
	cd ehd && make itags
	cd examples && make itags
	cd test && make itags
	make itags

etags:
	etags *.h grid/*.h

checklinks:
	$(LINKCHECKER) 	$(BASILISK_URL)/src/README 		\
			$(BASILISK_URL)/src/test/README 	\
			$(BASILISK_URL)/src/examples/README | 	\
		tee checklinks.log

checklinksfast:
	wget --spider -nd -nv -r \
		--reject-regex '.*[?]changes=.*' \
		--reject-regex '.*[?]history' \
		$(BASILISK_URL) 2>&1 | \
		grep -v ^unlink: | tee checklinks.log

changelog:
	darcs changes > ChangeLog

dist:
	darcs dist

diff:
	cd .. && tar czvf src/diff.tgz `darcs whatsnew -s | \
		sed 's/. .\/.*\/$$//g' | awk '{print $$2}'`

-include bview-server3D.d
-include bview-server2D.d

include $(BASILISK)/Makefile.defs
