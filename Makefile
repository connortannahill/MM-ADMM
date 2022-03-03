CC=g++
IDIRS=-I../src -I./lib/libigl/include/ -I./lib/eigen/ -I./lib/nanoflann/include/ -I./lib/json/single_include/nlohmann -I./lib/LASolver #-I./lib/LBFGSpp/include/
<<<<<<< HEAD
CFLAGS=-Wall -std=c++17 $(IDIRS) -O3 -msse2
# DEFS = -D THREADS
DEFS =
=======
# CFLAGS=-Wall -std=c++17 $(IDIRS) -O3 -msse2
CFLAGS=-Wall -std=c++17 $(IDIRS) -O3 -msse2 -fopenmp 
DEFS = -D THREADS
>>>>>>> a3f32770e569a4db99b8594c398d069d60676cac

SRC=main.cpp $(wildcard ./src/*.cpp) $(wildcard ./lib/LASolver/*.cpp)

mesh.exe : $(SRC) $(LIBS)
	$(CC) $(CFLAGS) $(DEFS) $^ -o $@

print-% : ;@echo $* = $($*)
