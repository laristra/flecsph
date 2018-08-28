#------------------------------------------------------------------------------#
# Copyright (c) 2018 Los Alamos National Security, LLC
# All rights reserved.
# ------------------------------------------------------------------------------#

# Makefile for generating static library from stellar_collapse

#CC = gcc
CC = mpicc
OUT_FILE_NAME1 = libTABEOS_SC.a
OUT_FILE_NAME2 = libTABEOS_RF.a

CFLAGS = -O2 -g -Wall -c 
LIBS = -lm -lmpi

OBJS1 = eos.o eos_gamma.o eos_poly.o util.o eos_stellar_collapse.o 
#OBJS1 = eos_gamma.o  
#OBJS2 = eos_gamma.o  
OBJS2 = eos.o eos_gamma.o eos_poly.o util.o root_finding.o

all: $(OUT_FILE_NAME1) $(OUT_FILE_NAME2)

$(OUT_FILE_NAME1): $(OBJS1)
	ar rcs $(OUT_FILE_NAME1) $(OBJS1)
	ranlib $(OUT_FILE_NAME1)

$(OUT_FILE_NAME2): $(OBJS2)
	ar rcs $(OUT_FILE_NAME2) $(OBJS2)
	ranlib $(OUT_FILE_NAME2)

eos.o : eos.c decs.h
	$(CC) $(CFLAGS) $(LIBS) eos.c
eos_gamma.o : eos_gamma.c decs.h
	$(CC) $(CFLAGS) $(LIBS) eos_gamma.c 
eos_poly.o : eos_poly.c decs.h
	$(CC) $(CFLAGS) $(LIBS) eos_poly.c 
eos_stella_collapse.o : eos_stellar_collapse.c eos_stella_collapse.h
	$(CC) $(CFLAGS) $(LIBS) eos_stellar_collapse.c 
root_finding.o : root_finding.c decs.h
	$(CC) $(CFLAGS) $(LIBS) root_finding.c 
util.o : util.c decs.h
	$(CC) $(CFLAGS) $(LIBS) util.c 

clean:
	rm -f *.o $(OUT_FILE_NAME1) $(OUT_FILE_NAME2)
