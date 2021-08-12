CC=g++
SRC=main.cpp $(wildcard ./src/*.cpp)
IDIRS=-I../src -I./lib/libigl/include/ -I./lib/eigen/ -I./lib/nanoflann/include/ #-I./lib/LBFGSpp/include/
CFLAGS=-Wall -std=c++11 $(IDIRS) -O3 -msse2 -pg
# LDFLAGS = -L/usr/local/opt/llvm/lib

SRC=main.cpp $(wildcard ./src/*.cpp)

mesh.exe : $(SRC) $(LIBS)
	$(CC) $(CFLAGS) $^ -o $@

print-% : ;@echo $* = $($*)
