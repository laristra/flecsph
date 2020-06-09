CXX=g++
CC=gcc
CXXFLAGS=-O3 --std=c++17
LDFLAGS=-lm -lhdf5

all: gen

.PHONY: clean
clean: 
	rm -rf nsID

gen: 
	${CXX} nsID_gen.cpp ${CXXFLAGS} -o nsID ${LDFLAGS}
