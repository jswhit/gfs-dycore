##intel fortran (jet,zeus)
#FC            = ifort
#OMPFLAG       = -openmp -no-ip
#FFLAGS        = -O3 -xHOST -assume byterecl -warn all -implicitnone -traceback  $(OMPFLAG)
#EXEC          = ../bin/gfs
#LDFLAGS       = $(OMPFLAG) -mkl=sequential
## jet
#LIBS          = -L/lfs1/projects/fim/whitaker/lib -lshtns -lfftw3
## zeus
#LIBS          = -L/scratch1/portfolios/BMC/fim/whitaker/lib -lshtns -lfftw3

#gnu fortran (mac)
FC = gfortran
OMPFLAG = -fopenmp
#FFLAGS=-O0 -pedantic -fimplicit-none -fbounds-check -fbacktrace -Wall -fcheck-array-temporaries -g
FFLAGS=-O3 -march=native -ffast-math -pipe -fomit-frame-pointer -fbacktrace $(OMPFLAG) 
EXEC          = ../bin/gfs
LDFLAGS       = $(OMPFLAG)  
LIBS          = -L/Users/jwhitaker/lib -lshtns -L/opt/local/lib -lfftw3 -llapack

.SUFFIXES: .o .f90 .F90

OBJS=	 params.o \
	 init.o\
	 kinds.o\
	 dyn_init.o\
	 dyn_run.o\
	 physcons.o\
	 sigio_module.o\
	 spectral_data.o\
	 pressure_data.o\
	 grid_data.o\
	 finalize.o\
	 dyn_finalize.o\
	 semimp_data.o\
	 run.o\
         phy_run.o\
         phy_init.o\
	 phy_finalize.o\
	 shtns.o


MAIN	= driver.f90

all: gfs

gfs: $(OBJS)
	$(FC) $(LDFLAGS) -o $(EXEC) $(MAIN) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) *.mod *genmod* $(EXEC)

.f90.o:
	$(FC) $(FFLAGS) -c $< 

.F90.o:
	$(FC) $(FFLAGS) -c $< 

.f.o:
	$(FC) $(FFLAGS) -free -c $< 

%.o : %.mod

kinds.o:	  kinds.f90

physcons.o:	  physcons.f90 kinds.o

constants.o:	  constants.f90 kinds.o

shtns.o:	  shtns.f90 kinds.o

sigio_module.o:   sigio_module.f90

params.o:	  params.f90 kinds.o sigio_module.o

init.o:  	  init.f90 params.o dyn_init.o phy_init.o kinds.o

finalize.o:  	  finalize.f90 params.o dyn_finalize.o phy_finalize.o

spectral_data.o:  spectral_data.f90 kinds.o

semimp_data.o:    semimp_data.f90 kinds.o shtns.o pressure_data.o physcons.o

grid_data.o:  	  grid_data.f90 kinds.o params.o

pressure_data.o:  pressure_data.f90 kinds.o params.o physcons.o

dyn_init.o:  	  dyn_init.f90 params.o sigio_module.o shtns.o spectral_data.o pressure_data.o physcons.o grid_data.o semimp_data.o

phy_init.o:       phy_init.f90 

phy_finalize.o:   phy_finalize.f90 

phy_run.o:        phy_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o

dyn_finalize.o:   dyn_finalize.f90 shtns.o spectral_data.o pressure_data.o grid_data.o

dyn_run.o:	  dyn_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o

run.o:   	  run.f90 dyn_run.o phy_run.o kinds.o params.o
