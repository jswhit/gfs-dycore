#intel fortran 
FC            = ifort
OMPFLAG       =  -openmp
FFLAGS        = -O3 -xHOST -assume byterecl -warn all -implicitnone -traceback  -convert big_endian $(OMPFLAG)
EXEC          = ../bin/gfs
EXEC_TEST     = ../bin/gfs_test
LDFLAGS       = $(OMPFLAG) -mkl=sequential 
# jet
LIBS          = -L/lfs1/projects/fim/whitaker/lib -lshtns -lfftw3 libgfs/libgfs.a w3lib/libw3.a
# zeus
#LIBS          = -L/scratch1/portfolios/BMC/fim/whitaker/lib -lshtns -lfftw3 libgfs/libgfs.a w3lib/libw3.a

##gnu fortran 
#FC = gfortran
#OMPFLAG = -fopenmp
##FFLAGS=-O0 -pedantic -fimplicit-none -fbounds-check -fbacktrace -Wall -fcheck-array-temporaries -g
#FFLAGS=-O3 -march=native -ffast-math -pipe -fomit-frame-pointer -fbacktrace $(OMPFLAG) 
#EXEC          = ../bin/gfs
#EXEC_TEST     = ../bin/gfs_test
#LDFLAGS       = $(OMPFLAG)  
#LIBS          = -L/Users/jwhitaker/lib -lshtns -lfftw3 -llapack libgfs/libgfs.a w3lib/libw3.a 
#LIBS_TEST     = -L/Users/jwhitaker/lib -lshtns -lfftw3 -llapack w3lib/libw3.a 

.SUFFIXES: .o .f90 .F90

OBJS=	 params.o \
	 init.o\
	 kinds.o\
	 dyn_init.o\
	 dyn_run.o\
	 physcons.o\
	 sigio_module.o\
	 sfcio_module.o\
	 phy_data.o\
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
	 funcphys.o\
	 patterngenerator.o\
	 stoch_data.o\
	 shtns.o


MAIN	= driver.f90

gfs: $(OBJS) 
	$(FC) $(LDFLAGS) -o $(EXEC) $(MAIN) $(OBJS) $(LIBS)

gfs-test: $(OBJS) 
	$(FC) $(LDFLAGS) -o $(EXEC_TEST) $(MAIN) $(OBJS) $(LIBS)

all: $(gfs)

libgfs : force_look
	cd libgfs; make clean; make

w3lib : force_look
	cd w3lib; make clean; make

force_look :
	true

clean:
	rm -f $(OBJS) *.mod *genmod* $(EXEC)
	ln -fs phy_init_gfs.f90 phy_init.f90
	ln -fs phy_data_gfs.f90 phy_data.f90
	ln -fs phy_finalize_gfs.f90 phy_finalize.f90
	ln -fs phy_run_gfs.f90 phy_run.f90

# Held-Suarez/Jablonowski-Williamson test
clean-test:
	rm -f $(OBJS) *.mod *genmod* $(EXEC_TEST)
	ln -fs phy_init_stub.f90 phy_init.f90
	ln -fs phy_data_stub.f90 phy_data.f90
	ln -fs phy_finalize_stub.f90 phy_finalize.f90
	ln -fs phy_run_hs.f90 phy_run.f90

.f90.o:
	$(FC) $(FFLAGS) -c $< 

.F90.o:
	$(FC) $(FFLAGS) -c $< 

.f.o:
	$(FC) $(FFLAGS) -free -c $< 

%.o : %.mod

kinds.o:	kinds.f90

physcons.o:	physcons.f90 kinds.o

constants.o:	constants.f90 kinds.o

shtns.o:	shtns.f90 kinds.o

sigio_module.o:  sigio_module.f90

sfcio_module.o:  sfcio_module.f90

params.o:	params.f90 kinds.o sigio_module.o

init.o:  	init.f90 params.o dyn_init.o phy_init.o phy_data.o physcons.o kinds.o 

finalize.o:  	finalize.f90 params.o dyn_finalize.o phy_finalize.o

spectral_data.o:  	spectral_data.f90 kinds.o

phy_data.o:  	phy_data.f90 kinds.o params.o sfcio_module.o physcons.o funcphys.o

semimp_data.o: semimp_data.f90 kinds.o shtns.o pressure_data.o physcons.o

grid_data.o:  	grid_data.f90 kinds.o params.o

pressure_data.o:  	pressure_data.f90 kinds.o params.o physcons.o

dyn_init.o:  	dyn_init.f90 params.o sigio_module.o shtns.o spectral_data.o pressure_data.o physcons.o grid_data.o semimp_data.o

phy_init.o:     phy_init.f90 phy_data.o

phy_finalize.o:     phy_finalize.f90 

phy_run.o:      phy_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o stoch_data.o
		$(FC) $(FFLAGS) -Ilibgfs -c phy_run.f90
	        
dyn_finalize.o:  dyn_finalize.f90 shtns.o spectral_data.o pressure_data.o grid_data.o

dyn_run.o:	dyn_run.f90 pressure_data.o shtns.o spectral_data.o physcons.o kinds.o params.o grid_data.o  stoch_data.o

run.o:   	run.f90 dyn_run.o phy_run.o kinds.o params.o stoch_data.o patterngenerator.o

pattergenerator.o:   	patterngenerator.f90 kinds.o shtns.o physcons.o

stoch_data.o  : stoch_data.f90 kinds.o params.o patterngenerator.o
