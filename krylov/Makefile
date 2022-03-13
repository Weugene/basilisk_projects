CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-L/usr/local/lib/ -L/home/e.sharaborin/opt/lib -L/trinity/home/e.sharaborin/opt/lib -L/trinity/home/e.sharaborin/opt/trinity/shared/opt/rh/gcc-8.3.0/lib64 -lcyaml -lm
C_INCLUDE_PATH=-I$(HOME)/opt/include
parallel=2

run:=mpirun -np $(parallel)  ./a.out > out 2> log
compile:= qcc -O2 -Wall -events -g  $(C_INCLUDE_PATH)
#compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -Wall -events -g  $(C_INCLUDE_PATH)
remove:=rm -rf .qcc*

all:
	$(remove)
	$(compile)  krylov.c $(CLIBS)
	./a.out
	#> out$(ind) 2> log$(ind)
