all:
	echo "choose target ring or ring_Popinet"
ring:
	$(eval CFILE=ring)
	$(eval ARGS=100 10 1 1 1 10)
	$(eval solver=a$(CFILE).out)
	$(eval dir=res_$(CFILE))
	$(eval FILENAME=rk_)
	@echo file:$(CFILE).c solver:$(solver) dir:$(dir) filename:$(FILENAME) args:$(ARGS)
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS) -o $(solver)
	mkdir -p $(dir)
	cp $(solver) $(dir)
	cd $(dir) && rm -rf ../.qcc* log out $(FILENAME)* snapshot-* nameYouWant*\ && ./$(solver) $(ARGS) > out 2> log

ring_Popinet_vof:
	$(eval CFILE=ring_Popinet)
	$(eval ARGS=100 10 1 1 1 9 1 )
	$(eval solver=a$(CFILE)_vof.out)
	$(eval dir=res_$(CFILE)_vof)
	$(eval FILENAME=rk_)
	@echo file:$(CFILE).c solver:$(solver) dir:$(dir) filename:$(FILENAME) args:$(ARGS)
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS) -o $(solver)
	mkdir -p $(dir)
	cp $(solver) $(dir)
	cd $(dir) && rm -rf ../.qcc* log out $(FILENAME)* snapshot-* nameYouWant*\ && ./$(solver) $(ARGS) > out 2> log

ring_Popinet_moving_cyl:
	$(eval CFILE=ring_Popinet)
	$(eval ARGS=100 10 1 1 1 9 0 )
	$(eval solver=a$(CFILE)_moving_cyl.out)
	$(eval dir=res_$(CFILE)_moving_cyl)
	$(eval FILENAME=rk_)
	@echo file:$(CFILE).c solver:$(solver) dir:$(dir) filename:$(FILENAME) args:$(ARGS)
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS) -o $(solver)
	mkdir -p $(dir)
	cp $(solver) $(dir)
	cd $(dir) && rm -rf ../.qcc* log out $(FILENAME)* snapshot-* nameYouWant*\ && ./$(solver) $(ARGS) > out 2> log
