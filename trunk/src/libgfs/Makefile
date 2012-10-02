include ../Makefile.conf
EXEC	 = ./libgfs.a

OPTSB=$(FFLAGS) $(AUTOFLAG)

OPTSM=$(OPTSB)
OPTSIOM=$(OPTSB)
OPTS_SERM=$(OPTSB) $(ARCHM)
OPTS90M=$(OPTSB)  $(FREEFLAG)
OPTS90AM=$(OPTSB) $(FREEFLAG)
 
OPTS     = $(OPTSM)
OPTS_SER = $(OPTS_SERM)
OPTS90   = $(OPTS90M)
OPTS90A  = $(OPTS90AM)
OPTSIO	 = $(OPTSIOM)
FFLAG90  = $(OPTS90) 
FFLAG90A = $(OPTS90A) 
FFLAGX   = $(OPTS) $(TRAPS) 
FFLAGIO  = $(OPTSIO) $(TRAPS) 
FFLAGY   = $(OPTS) 
FFLAGM   = $(OPTS) $(TRAPS) $(DEBUG)  
FFLAGM3  = $(OPTS) $(TRAPS) $(DEBUG) $(FREEFLAG)
FFLAGSF  = $(OPTSB) $(FREEFLAG) 
FFLAGSI  = $(OPTSB) $(FREEFLAG)
FFLAGB   = $(OPTSB) 
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
          phy_data_gfs.o \
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
	$(FC) $(FFLAGX) -c -d $< 
.f.o:
	$(FC) $(FFLAGX) -c $< 

# *****************************************************************
#
cnvc90.o:	cnvc90.f
		$(FC) $(FFLAGM) -c cnvc90.f

calpreciptype.o:	calpreciptype.f
		$(FC) $(FFLAGM3) -c calpreciptype.f

dcyc2.o:	dcyc2.f
		$(FC) $(FFLAGM) -c dcyc2.f

dcyc2.pre.rad.o:	dcyc2.pre.rad.f
		$(FC) $(FFLAGM) -c dcyc2.pre.rad.f

gbphys.o:	gbphys.f
		$(FC) $(FFLAGM) -c gbphys.f

get_prs.o:	get_prs.f
		$(FC) $(FFLAGM) -c get_prs.f

gscond.o:	gscond.f
		$(FC) $(FFLAGM) -c gscond.f

gsmddrive.o:	gsmddrive.f
		$(FC) $(FFLAGM) -c gsmddrive.f

module_bfmicrophysics.o:	module_bfmicrophysics.f
		$(FC) $(FFLAGM) -c module_bfmicrophysics.f

gwdps.o:	gwdps.f
		$(FC) $(FFLAGM) -c gwdps.f

gwdc.o:		gwdc.f
		$(FC) $(FFLAGM) -c gwdc.f

hyb2press.o:	hyb2press.f
		$(FC) $(FFLAGM) -c hyb2press.f

lrgsclr.o:	lrgsclr.f
		$(FC) $(FFLAGM) -c lrgsclr.f

moninp.o:	moninp.f
		$(FC) $(FFLAGM) -c moninp.f

moninp1.o:	moninp1.f
		$(FC) $(FFLAGM) -c moninp1.f

moninq.o:	moninq.f
		$(FC) $(FFLAGM) -c moninq.f

moninq1.o:	moninq1.f
		$(FC) $(FFLAGM) -c moninq1.f

mstadb.o:	mstadb.f
		$(FC) $(FFLAGM) -c mstadb.f

mstcnv.o:	mstcnv.f
		$(FC) $(FFLAGM) -c mstcnv.f

ozphys.o:	ozphys.f
		$(FC) $(FFLAGM) -c ozphys.f

precpd.o:	precpd.f
		$(FC) $(FFLAGM) -c precpd.f

rascnvv2.o:	rascnvv2.f
		$(FC) $(FFLAGM) -c rascnvv2.f

sfc_sice.o:	sfc_sice.f
		$(FC) $(FFLAGM) -c sfc_sice.f

sfc_land.o:	sfc_land.f
		$(FC) $(FFLAGM) -c sfc_land.f
                                                                            
progt2.o:	progt2.f
		$(FC) $(FFLAGM) -c progt2.f

sfc_drv.o:	sfc_drv.f
		$(FC) $(FFLAGM) -c sfc_drv.f

sflx.o: 	sflx.f phy_data_gfs.o
		$(FC) $(FFLAGM) -c sflx.f

sfc_ocean.o:	sfc_ocean.f
		$(FC) $(FFLAGM) -c sfc_ocean.f

sfc_nst.o:	sfc_nst.f
		$(FC) $(FFLAGM) -c sfc_nst.f

sfc_diff.o:	sfc_diff.f
		$(FC) $(FFLAGM) -c sfc_diff.f

sfc_diag.o:	sfc_diag.f
		$(FC) $(FFLAGM) -c sfc_diag.f

sascnv.o:	sascnv.f
		$(FC) $(FFLAGM) -c sascnv.f

sascnvn.o:	sascnvn.f
		$(FC) $(FFLAGM) -c sascnvn.f

tridi2t3.o:	tridi2t3.f
		$(FC) $(FFLAGM) -c tridi2t3.f

shalcv.o:	shalcv.f
		$(FC) $(FFLAGM) -c shalcv.f

shalcv_opr.o:	shalcv_opr.f
		$(FC) $(FFLAGM) -c shalcv_opr.f

shalcnv.o:	shalcnv.f
		$(FC) $(FFLAGM) -c shalcnv.f

# *****************************************************************

radlw_param.o:	radlw_param.f
		$(FC) $(FFLAGM) -c radlw_param.f


radlw_datatb.o:	radlw_datatb.f
		$(FC) $(FFLAGM) -c radlw_datatb.f


radlw_main.o:	radlw_main.f
		$(FC) $(FFLAGM) -c radlw_main.f

radsw_param.o:	radsw_param.f
		$(FC) $(FFLAGM) -c radsw_param.f

radsw_datatb.o:	radsw_datatb.f
		$(FC) $(FFLAGM) -c radsw_datatb.f

radsw_main.o:	radsw_main.f
		$(FC) $(FFLAGM) -c radsw_main.f

radiation_astronomy.o:	radiation_astronomy.f
		$(FC) $(FFLAGM) -c radiation_astronomy.f

radiation_aerosols.o:	radiation_aerosols.f
		$(FC) $(FFLAGM) -c radiation_aerosols.f

radiation_gases.o:	radiation_gases.f
		$(FC) $(FFLAGM) -c radiation_gases.f

radiation_clouds.o:	radiation_clouds.f
		$(FC) $(FFLAGM) -c radiation_clouds.f

radiation_surface.o:	radiation_surface.f
		$(FC) $(FFLAGM) -c radiation_surface.f

grrad.o:	grrad.f
		$(FC) $(FFLAGM) -c grrad.f

progtm_module.o:	progtm_module.f
		$(FC) $(FFLAGM) -c progtm_module.f

machine.o:	machine.f kinds.o
		$(FC) $(FFLAGM)  -c machine.f

mpi_quit.o:	mpi_quit.f90
		$(FC) $(FFLAG90)  -c mpi_quit.f90

kinds.o:	../kinds.f90
		$(FC) $(FFLAG90) -c ../kinds.f90

#
# *****************************************************************
#
physcons.o:	../physcons.f90
		$(FC) $(FFLAG90) -c ../physcons.f90

phy_data_gfs.o:	../phy_data_gfs.f90
		$(FC) $(FFLAG90) -c ../phy_data_gfs.f90

params.o:	../params.f90
		$(FC) $(FFLAG90) -c ../params.f90

funcphys.o:	../funcphys.f90
		$(FC) $(FFLAG90) -c ../funcphys.f90

sigio_module.o:	../sigio_module.f90
		$(FC) $(FFLAG90) -c ../sigio_module.f90

sfcio_module.o:	../sfcio_module.f90
		$(FC) $(FFLAG90) -c ../sfcio_module.f90

shtns.o:	../shtns.f90
		$(FC) $(FFLAG90) -c ../shtns.f90

nstio_module.o:	 nstio_module.f
		$(FC) $(FFLAG90) -c nstio_module.f

module_nst_parameters.o:	module_nst_parameters.f
		$(FC) $(FFLAG90) -c module_nst_parameters.f

module_nst_water_prop.o:	 module_nst_water_prop.f
		$(FC) $(FFLAG90) -c module_nst_water_prop.f

module_nst_model.o:		module_nst_model.f
		$(FC) $(FFLAG90) -c module_nst_model.f

gribit.o:	gribit.f
		$(FC) $(FFLAGIO) -c gribit.f
