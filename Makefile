CC=g++
SRC=main.cpp $(wildcard ./src/*.cpp)
IDIRS=-I../src -I./lib/eigen/ -I./lib/nanoflann/include/
CFLAGS=-Wall -std=c++11 $(IDIRS) -O3 -msse2 -fopenmp -D THREADS
# CFLAGS=-Wall -std=c++11 $(IDIRS) -O3 -msse2 -fopenmp
# DEFS = -D THREADS
DEFS =

SRC=main.cpp $(wildcard ./src/*.cpp)

mesh.exe : $(SRC) $(LIBS)
	$(CC) $(CFLAGS) $(DEFS) $^ -o $@

print-% : ;@echo $* = $($*)
