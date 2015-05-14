# Makefile for the program mpi_ex_send.c
objects = mpi_ex_send.o
EXE = mpi_ex_send
CFLAGS ?= -O3

mpi_ex_send: $(objects)
	cc -o mpi_ex_send $(objects) -lm

mpi_ex_send.o: mpi_ex_send.c mpi_ex_send.h

all: mpi_ex_send

clean:
	rm $(EXE) $(objects)

