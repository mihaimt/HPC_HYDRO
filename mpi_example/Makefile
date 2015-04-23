# Makefile for the program mpi_example.c
objects = mpi_example.o
EXE = mpi_example
CFLAGS ?= -O3

mpi_example: $(objects)
	cc -o mpi_example $(objects) -lm

mpi_example.o: mpi_example.c mpi_example.h

all: mpi_example

clean:
	rm $(EXE) $(objects)

