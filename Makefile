CC=g++
IDIRS=-I../src -I./lib/libigl/include/ -I./lib/eigen/ -I./lib/nanoflann/include/ -I./lib/json/single_include/nlohmann -I./lib/LASolver #-I./lib/LBFGSpp/include/
CFLAGS=-Wall -std=c++17 $(IDIRS) -O3 -msse2
# DEFS = -D THREADS
DEFS =

SRC=main.cpp $(wildcard ./src/*.cpp) $(wildcard ./lib/LASolver/*.cpp)

mesh.exe : $(SRC) $(LIBS)
	$(CC) $(CFLAGS) $(DEFS) $^ -o $@

print-% : ;@echo $* = $($*)
