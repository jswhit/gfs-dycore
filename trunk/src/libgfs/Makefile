EXEC	 = ./libgfs.a

F77=ifort
F90=ifort
FRRM=-FR
OPTSB=-O3 -convert big_endian -xHOST -recursive
OPTSBT=$(OPTSB) -traceback

#F77=gfortran
#F90=gfortran
#FRRM=-ffree-form
#OPTSB=-O3 -march=native -frecursive -ffast-math -pipe -fomit-frame-pointer -fno-range-check
#OPTSBT=$(OPTSB) -fbacktrace

OPTSM=$(OPTSBT)
OPTSIOM=$(OPTSBT)
OPTS_SERM=$(OPTSBT) $(ARCHM)
OPTS90M=$(OPTSBT)  $(FRRM)
OPTS90AM=$(OPTSBT) $(FRRM)
 
OPTS     = $(OPTSM)
OPTS_SER = $(OPTS_SERM)
OPTS90   = $(OPTS90M)
OPTS90A  = $(OPTS90AM)
OPTSIO	 = $(OPTSIOM)
FFLAG90  = $(OPTS90) 
FFLAG90A = $(OPTS90A) 
FFLAGS   = $(OPTS) $(TRAPS)  
FFLAGX   = $(OPTS) $(TRAPS) 
FFLAGIO  = $(OPTSIO) $(TRAPS) 
FFLAGY   = $(OPTS) 
FFLAGM   = $(OPTS) $(TRAPS) $(DEBUG)  
FFLAGM3  = $(OPTS) $(TRAPS) $(DEBUG) $(FRRM)
FFLAGSF  = $(OPTSBT) $(FRRM) 
FFLAGSI  = $(OPTSBT) $(FRRM)
FFLAGB   = $(OPTSBT) 
FFLAG_SER = $(OPTS_SER) $(DOUBLEFLAG) $(TRAPS) $(DEBUG) 
FFLAGSC   =

LDR	 = $(LDRM)
LDFLAGS = 
LIBS =

.SUFFIXES: .o .f .F .h


.SUFFIXES: .o .f .F .h
#
# *****************************************************************
#

OBJ_MOD	= kinds.o             \
	  mpi_quit.o	\
	  machine.o \
          physcons.o \
          sigio_module.o \
          sfcio_module.o \
          shtns.o \
          params.o \
          phy_data.o \
	  funcphys.o \
	  iounitdef.o           \
	  progtm_module.o       \
	  rascnvv2.o            \
	  mersenne_twister.o    \
	  tracer_const_h-new.o \
	  module_bfmicrophysics.o  


OBJS_RAD	=              \
	radlw_param.o          \
	radlw_datatb.o         \
	radlw_main.o           \
	radsw_param.o          \
	radsw_datatb.o         \
	radsw_main.o           \
	radiation_astronomy.o  \
	radiation_aerosols.o   \
	radiation_gases.o      \
	radiation_clouds.o     \
	radiation_surface.o    \
	grrad.o                

OBJS_PHY= \
	ozphys.o            \
	calpreciptype.o     \
	gbphys.o            \
	dcyc2.o             \
	dcyc2.pre.rad.o     \
	sfc_drv.o           \
	sfc_land.o          \
	progt2.o            \
	sfc_sice.o          \
	sfc_ocean.o         \
	sfc_diff.o          \
	sfc_diag.o          \
	sflx.o              \
	moninp.o            \
	moninp1.o           \
	moninq.o            \
	moninq1.o           \
	tridi2t3.o          \
	gwdps.o             \
	gwdc.o              \
	sascnv.o            \
	sascnvn.o           \
	cnvc90.o            \
	shalcv.o            \
	shalcv_opr.o        \
	shalcnv.o           \
	lrgsclr.o           \
	gscond.o            \
	precpd.o            \
	mstadb.o            \
	mstadbtn.o          \
	mstcnv.o            \
	get_prs.o           \
        gribit.o	\
	gsmddrive.o         

SRC	= $(OBJ_MOD:.o=.f) $(OBJS_RAD:.o=.f) $(OBJS_PHY:.o=.f)

#
# *****************************************************************
#
all: model-mpi

model-mpi: $(OBJ_MOD) $(OBJS_PHY) $(OBJS_RAD) 
	ar -ruv $(EXEC) $(OBJ_MOD) $(OBJS_PHY) $(OBJS_RAD) 

clean:
	rm -f *.o *.mod $(EXEC)

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 

# *****************************************************************
#
cnvc90.o:	cnvc90.f
		$(F77) $(FFLAGM) -c cnvc90.f

calpreciptype.o:	calpreciptype.f
		$(F77) $(FFLAGM3) -c calpreciptype.f

dcyc2.o:	dcyc2.f
		$(F77) $(FFLAGM) -c dcyc2.f

dcyc2.pre.rad.o:	dcyc2.pre.rad.f
		$(F77) $(FFLAGM) -c dcyc2.pre.rad.f

gbphys.o:	gbphys.f
		$(F77) $(FFLAGM) -c gbphys.f

get_prs.o:	get_prs.f
		$(F77) $(FFLAGM) -c get_prs.f

gscond.o:	gscond.f
		$(F77) $(FFLAGM) -c gscond.f

gsmddrive.o:	gsmddrive.f
		$(F77) $(FFLAGM) -c gsmddrive.f

module_bfmicrophysics.o:	module_bfmicrophysics.f
		$(F77) $(FFLAGM) -c module_bfmicrophysics.f

gwdps.o:	gwdps.f
		$(F77) $(FFLAGM) -c gwdps.f

gwdc.o:		gwdc.f
		$(F77) $(FFLAGM) -c gwdc.f

hyb2press.o:	hyb2press.f
		$(F77) $(FFLAGM) -c hyb2press.f

lrgsclr.o:	lrgsclr.f
		$(F77) $(FFLAGM) -c lrgsclr.f

moninp.o:	moninp.f
		$(F77) $(FFLAGM) -c moninp.f

moninp1.o:	moninp1.f
		$(F77) $(FFLAGM) -c moninp1.f

moninq.o:	moninq.f
		$(F77) $(FFLAGM) -c moninq.f

moninq1.o:	moninq1.f
		$(F77) $(FFLAGM) -c moninq1.f

mstadb.o:	mstadb.f
		$(F77) $(FFLAGM) -c mstadb.f

mstcnv.o:	mstcnv.f
		$(F77) $(FFLAGM) -c mstcnv.f

ozphys.o:	ozphys.f
		$(F77) $(FFLAGM) -c ozphys.f

precpd.o:	precpd.f
		$(F77) $(FFLAGM) -c precpd.f

rascnvv2.o:	rascnvv2.f
		$(F77) $(FFLAGM) -c rascnvv2.f

sfc_sice.o:	sfc_sice.f
		$(F77) $(FFLAGM) -c sfc_sice.f

sfc_land.o:	sfc_land.f
		$(F77) $(FFLAGM) -c sfc_land.f
                                                                            
progt2.o:	progt2.f
		$(F77) $(FFLAGM) -c progt2.f

sfc_drv.o:	sfc_drv.f
		$(F77) $(FFLAGM) -c sfc_drv.f

sflx.o: 	sflx.f phy_data.o
		$(F77) $(FFLAGM) -c sflx.f

sfc_ocean.o:	sfc_ocean.f
		$(F77) $(FFLAGM) -c sfc_ocean.f

sfc_nst.o:	sfc_nst.f
		$(F77) $(FFLAGM) -c sfc_nst.f

sfc_diff.o:	sfc_diff.f
		$(F77) $(FFLAGM) -c sfc_diff.f

sfc_diag.o:	sfc_diag.f
		$(F77) $(FFLAGM) -c sfc_diag.f

sascnv.o:	sascnv.f
		$(F77) $(FFLAGM) -c sascnv.f

sascnvn.o:	sascnvn.f
		$(F77) $(FFLAGM) -c sascnvn.f

tridi2t3.o:	tridi2t3.f
		$(F77) $(FFLAGM) -c tridi2t3.f

shalcv.o:	shalcv.f
		$(F77) $(FFLAGM) -c shalcv.f

shalcv_opr.o:	shalcv_opr.f
		$(F77) $(FFLAGM) -c shalcv_opr.f

shalcnv.o:	shalcnv.f
		$(F77) $(FFLAGM) -c shalcnv.f

# *****************************************************************

radlw_param.o:	radlw_param.f
		$(F77) $(FFLAGM) -c radlw_param.f


radlw_datatb.o:	radlw_datatb.f
		$(F77) $(FFLAGM) -c radlw_datatb.f


radlw_main.o:	radlw_main.f
		$(F77) $(FFLAGM) -c radlw_main.f

radsw_param.o:	radsw_param.f
		$(F77) $(FFLAGM) -c radsw_param.f

radsw_datatb.o:	radsw_datatb.f
		$(F77) $(FFLAGM) -c radsw_datatb.f

radsw_main.o:	radsw_main.f
		$(F77) $(FFLAGM) -c radsw_main.f

radiation_astronomy.o:	radiation_astronomy.f
		$(F77) $(FFLAGM) -c radiation_astronomy.f

radiation_aerosols.o:	radiation_aerosols.f
		$(F77) $(FFLAGM) -c radiation_aerosols.f

radiation_gases.o:	radiation_gases.f
		$(F77) $(FFLAGM) -c radiation_gases.f

radiation_clouds.o:	radiation_clouds.f
		$(F77) $(FFLAGM) -c radiation_clouds.f

radiation_surface.o:	radiation_surface.f
		$(F77) $(FFLAGM) -c radiation_surface.f

grrad.o:	grrad.f
		$(F77) $(FFLAGM) -c grrad.f

progtm_module.o:	progtm_module.f
		$(F77) $(FFLAGM) -c progtm_module.f

machine.o:	machine.f kinds.o
		$(F77) $(FFLAGM)  -c machine.f

mpi_quit.o:	mpi_quit.f90
		$(F90) $(FFLAG90)  -c mpi_quit.f90

kinds.o:	../kinds.f90
		$(F90) $(FFLAG90) -c ../kinds.f90

#
# *****************************************************************
#
physcons.o:	../physcons.f90
		$(F90) $(FFLAG90) -c ../physcons.f90

phy_data.o:	../phy_data.f90
		$(F90) $(FFLAG90) -c ../phy_data.f90

params.o:	../params.f90
		$(F90) $(FFLAG90) -c ../params.f90

funcphys.o:	../funcphys.f90
		$(F90) $(FFLAG90) -c ../funcphys.f90

sigio_module.o:	../sigio_module.f90
		$(F77) $(FFLAG90) -c ../sigio_module.f90

sfcio_module.o:	../sfcio_module.f90
		$(F77) $(FFLAG90) -c ../sfcio_module.f90

shtns.o:	../shtns.f90
		$(F77) $(FFLAG90) -c ../shtns.f90

nstio_module.o:	 nstio_module.f
		$(F77) $(FFLAG90) -c nstio_module.f

module_nst_parameters.o:	module_nst_parameters.f
		$(F77) $(FFLAG90) -c module_nst_parameters.f

module_nst_water_prop.o:	 module_nst_water_prop.f
		$(F77) $(FFLAG90) -c module_nst_water_prop.f

module_nst_model.o:		module_nst_model.f
		$(F77) $(FFLAG90) -c module_nst_model.f

gribit.o:	gribit.f
		$(F77) $(FFLAGIO) -c gribit.f
