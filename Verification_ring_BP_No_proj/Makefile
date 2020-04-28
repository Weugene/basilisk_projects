all:
	echo "choose target ring or ring_Popinet"
ring:
	$(eval CFILE=ring)
	$(eval ARGS=100 0.01 1 1 1 10)
	$(eval solver=a$(CFILE).out)
	$(eval dir=res_$(CFILE))
	$(eval FILENAME=rk_)
	@echo file:$(CFILE).c solver:$(solver) dir:$(dir) filename:$(FILENAME) args:$(ARGS)
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS) -o $(solver)
	mkdir -p $(dir)
	cp $(solver) $(dir)
	cd $(dir) && rm -rf ../.qcc* log out $(FILENAME)* snapshot-* nameYouWant*\ && ./$(solver) $(ARGS)> out 2> log

ring_Popinet:
	$(eval CFILE=ring_Popinet)
	$(eval ARGS=100 0.01 1 1 1 10 1)
	$(eval solver=a$(CFILE).out)
	$(eval dir=res_$(CFILE)_dir)
	$(eval FILENAME=rk_)
	@echo file:$(CFILE).c solver:$(solver) dir:$(dir) filename:$(FILENAME) args:$(ARGS)
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS) -o $(solver)
	mkdir -p $(dir)
	cp $(solver) $(dir)
	cd $(dir) && rm -rf ../.qcc* log out $(FILENAME)* snapshot-* nameYouWant*\ && ./$(solver) $(ARGS)> out 2> log
