# Makefile for the program mpi_ex_isend_irecv.c
objects = mpi_ex_isend_irecv.o
EXE = mpi_ex_isend_irecv
CFLAGS ?= -O3 

mpi_ex_isend_irecv: $(objects)
	cc -o mpi_ex_isend_irecv $(objects) -lm

mpi_ex_isend_irecv.o: mpi_ex_isend_irecv.c mpi_ex_isend_irecv.h

all: mpi_ex_isend_irecv

clean:
	rm $(EXE) $(objects)

