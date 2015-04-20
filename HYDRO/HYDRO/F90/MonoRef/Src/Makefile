# Link make.inc with the targeted architecture in Arch directory
# For example : ln -s make_sp6 make.inc
#--------------------------------------
include ../../../Arch/make.inc

# Sources directories
#--------------------
SRCDIR = .

# Objects
#--------
OBJ = module_hydro_commun.o module_hydro_IO.o module_hydro_utils.o module_hydro_principal.o main.o

default: mono

mono: FCOMP = $(F90)
mono: FCOMPFLAGS = $(FFLAGS)
mono: LINKFLAGS = $(LDFLAGS)
mono: $(OBJ)
	$(FCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mono

mpi: FCOMP = $(MPIF90)
mpi: FCOMPFLAGS = $(FFLAGS) $(FFLAGS_MPI)
mpi: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_MPI)
mpi: $(OBJ)
	$(FCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mpi

omp: FCOMP = $(F90)
omp: FCOMPFLAGS = $(FFLAGS) $(FFLAGS_OMP)
omp: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_OMP)
omp: $(OBJ)
	$(FCOMP) $(LINKFLAGS) $(OBJ) -o hydro_omp

mpiomp: FCOMP = $(MPIF90)
mpiomp: FCOMPFLAGS = $(FFLAGS) $(FFLAGS_MPI) $(FFLAGS_OMP)
mpiomp: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_MPI) $(LDFLAGS_OMP)
mpiomp: $(OBJ)
	$(FCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mpiomp

clean:
	rm -f *.o *.mod core* 

cleanall:
	rm -f *.o *.mod core* hydro_mono hydro_mpi hydro_omp hydro_mpiomp

%.o:    $(SRCDIR)/%.f90
	$(FCOMP) $(FCOMPFLAGS) -c $<

# Dependencies
#-------------
main.o : main.f90 module_hydro_principal.o module_hydro_IO.o module_hydro_commun.o

module_hydro_commun.o : module_hydro_commun.f90

module_hydro_IO.o : module_hydro_IO.f90 module_hydro_commun.o

module_hydro_principal.o : module_hydro_principal.f90 module_hydro_commun.o module_hydro_utils.o

module_hydro_utils.o : module_hydro_utils.f90 module_hydro_commun.o
