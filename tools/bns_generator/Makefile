
all: gen

clean: 
	rm -rf bns_generator

gen: ; mpiicpc generator_bns_3D.cpp -O3 -o bns_generator -lm -I ~/local/hdf5_intel/include/ -L ~/local/hdf5_intel/lib/ --std=c++17 -lhdf5
