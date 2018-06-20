CXX=g++
CC=gcc
CXXFLAGS=-O3 --std=c++17
LDFLAGS=-lm -lhdf5

all: gen

.PHONY: clean
clean: 
	rm -rf bns_generator

gen: 
	${CXX} generator_bns_3D.cpp ${CXXFLAGS} -o bns_generator ${LDFLAGS}
