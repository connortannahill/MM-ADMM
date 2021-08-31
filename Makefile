CC=g++
SRC=main.cpp $(wildcard ./src/*.cpp)
IDIRS=-I../src -I./lib/libigl/include/ -I./lib/eigen/ -I./lib/nanoflann/include/ #-I./lib/LBFGSpp/include/
CFLAGS=-Wall -std=c++11 $(IDIRS) -O3 -msse2 -fopenmp
DEFS = -D NUMTHREADS=2

SRC=main.cpp $(wildcard ./src/*.cpp)

mesh.exe : $(SRC) $(LIBS)
	$(CC) $(CFLAGS) $(DEFS) $^ -o $@

# Parallel test compilation
# meshP1.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=1 $^ -o $@

# meshP2.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=2 $^ -o $@

# meshP4.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=4 $^ -o $@

# meshP8.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=8 $^ -o $@

# meshP16.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=16 $^ -o $@

# meshP32.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=32 $^ -o $@

# meshP64.exe : $(SRC) $(LIBS)
# 	$(CC) $(CFLAGS) -D NUMTHREADS=64 $^ -o $@

# par = meshP1.exe meshP2.exe meshP4.exe meshP8.exe meshP16.exe meshP32.exe meshP64.exe
# par: $(par)


print-% : ;@echo $* = $($*)
