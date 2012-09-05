module params
! holds model parameters
! Public subroutines:
! read_namelist: read namelist
 use kinds, only: r_kind,r_double
 use sigio_module, only: sigio_head, sigio_srhead, sigio_sropen, sigio_sclose

 implicit none
 private

 public :: read_namelist,initfile,sfcinitfile,fhmax,dt,ntmax,ndimspec,nlons,nlats,&
 tstart,ndiss,efold,nlevs,ntrunc,sighead,dry,explicit,heldsuarez,jablowill,&
 ntout,fhdfi,fhout,fhzer,idate_start,adiabatic,hdif_fac,hdif_fac2,fshk,ntrac,ntoz,ntclw,&
 postphys,timestepsperhr,ncw,taustratdamp,polar_opt,ntdfi,&
! gfs phys parameters.
 nmtvr,fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,&
 iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,psautco,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,mstrat,&
 pre_rad,bkgd_vdif_m,bkgd_vdif_s,bkgd_vdif_h,gloopb_filter

 character(len=500) :: initfile ! init cond filename
 character(len=500) :: sfcinitfile ! surface init cond filename
 integer            :: fhmax ! hours to run
 integer            :: fhout ! interval for IO
 integer            :: fhzer ! interval to zero accumulated arrays
! half window length (hrs) for digital filter launch (=0 mean no dfi)
 integer            :: fhdfi=0
 real(r_double)     :: deltim=0    ! namelist input time step (secs)
 real(r_double)     :: dt    ! time step (secs) (=deltim or 3600/timestepsperhr)
 integer    :: ntmax ! time steps to run
 integer    :: ntdfi ! number of time steps in dfi window is 2*ntdfi+1
 integer    :: nlons ! number of longitudes on grid
 integer    :: nlats ! number of latitudes on grid
 integer    :: nlevs ! number of levels on grid
 integer    :: ntrunc ! spectral truncation
 integer    :: ndimspec ! spectral array dimension
 type(sigio_head),save  :: sighead ! header struct from initfile
 logical    :: dry = .false. ! no moisture, cloud condensate or ozone.
 logical    :: adiabatic = .false. ! don't call physics
 ! held-suarez forcing
 logical    :: heldsuarez = .false.
 ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
 ! idealized baroclinic instability test case.
 logical    :: jablowill = .false.
 ! use explicit time differencing
 ! if .true., explicit RK3 is used.
 ! if .false., semi-implicit RK3 (Kar, 2006, MWR p. 2916,
 ! http://dx.doi.org/10.1175/MWR3214.1) is used.
 ! explicit or semi-implicit rk3
 !logical    :: explicit = .true. ! use explicit rk3
 logical    :: explicit = .false. ! use semi-implicit rk3
 ! starting forecast time in seconds (read in from initfile)
 real(r_kind) :: tstart
 integer    :: idate_start(4) ! starting date (hr,month,day,year)
 integer    :: ntout ! time step interval for IO
 integer    :: ndiss=0 ! hyperdiffusion order (0 for GFS defaults)
 real(r_double) :: polar_opt=1.e-10 ! polar optimization threshold for transforms
 ! efolding scale for smallest resolvable wave (0 for GFS defaults)
 real(r_kind) :: efold=0. 
 real(r_kind) :: hdif_fac=1.0 ! multiplier for height-dep part of hyper-diff
 real(r_kind) :: hdif_fac2=1.0 ! multiplier to increase hyper-diffusion
 ! amplitude of vertically varying part of hyper-diff (1 means no variation,
 ! zero gives GFS resolution dependent defaults)
 real(r_kind) :: fshk=0 
 integer :: ntrac=3 ! number of tracers (including specific humidity)
 integer :: ntoz=2 ! ozone tracer number
 integer :: ntclw=3 ! cloud condensate tracer number
 integer :: nmtvr=14 ! number of fields in mtnvar file.
 real(r_kind) :: taustratdamp=5.*86400. ! extra linear drag near top of model
! parameters relevant for GFS physics
 ! interval in hours to call long-wave radiation (0 means every time step)
 real(r_kind) :: fhlwr=0 
 ! interval in hours to call short-wave radiation (0 means every time step)
 real(r_kind) :: fhswr=0 
 ! ictm controls source for controls source for time sensitive external data (e.g. CO2,
 ! solcon, aerosols, etc)
 integer :: ictm=0 ! use data at initial cond time, or latest available.
 ! isol controls solar constant data source.
 integer :: isol=0 ! use prescribed solar constant
 ! ico2 controls co2 data source for radiation
 integer :: ico2=0 ! prescribed global mean co2
 ! iaer controls aerosols scheme selection
 integer :: iaer=1 ! default aerosol
 ! ialb controls surface albedo for sw radiation
 integer :: ialb=0 ! use climo albedo based on sfc type
 ! iems controls surface emissivity and sfc air/ground temp for lw radiation
 integer :: iems=0 ! used fixed value of 1.0
 ! isubc_sw,lw controls sub-column cloud approximation in radiation
 integer :: isubc_sw=0 ! sw clouds without sub-column approximation
 integer :: isubc_lw=0 ! lw clouds without sub-column approximation
 ! iovr_sw,lw controls cloud overlapping method in radiation
 integer :: iovr_sw=1 ! max-random overlap clouds
 integer :: iovr_lw=1 ! max-random overlap clouds
 logical :: newsas=.true. ! 'new. SAS convection scheme
 logical :: ras=.false. ! RAS convection scheme.
 logical :: sashal=.true. ! 'new' mass-sflux based shallow convection
 logical :: mstrat=.false. ! flag for moorthi approach for stratus
 logical :: pre_rad=.false. ! flag for debugging in gbphys
 integer :: num_p3d=4 ! # of 3d microphysics fields (value for Zhao microphysics)
 integer :: num_p2d=3 ! # of 2d microphysics fields (value for Zhao)
 logical :: crick_proof = .false.
 logical :: ccnorm = .false.
 logical :: norad_precip = .false. ! only used for Ferrier microphysics
 ! range of droplet number concentrations for Ferrier scheme.
 integer :: ncw(2)=(/50,150/)
 ! critical relative humidity at the surface, PBL and atmosphere top   
 real(r_kind) :: crtrh(3)=(/0.85,0.85,0.85/)
 ! multiplication factors for cdmb and gwd  (Mtn Blking and GWD tuning factors)
 real(r_kind) :: cdmbgwd(2)=(/1.0,1.0/)
 ! RAS convection mult factor for critical cloud work function.
 real(r_kind) :: ccwf(2)=(/1.0,1.0/)
 ! factor for cloud condensate detrainment from cloud edges in RAS.
 real(r_kind) :: dlqf(2)=(/0.0,0.0/)
 ! critical cloud top entrainment instability criteria (for mstrat=.true.)
 real(r_kind) :: ctei_rm(2)=(/10.0,10.0/)
 ! auto conversion coeff from ice to snow for Zhao microphysics
 real(r_kind) :: psautco(2)=(/4.0e-4,4.0e-4/)
 ! auto conversion coeff from cloud to rain for Zhao microphysics
 real(r_kind) :: prautco(2)=(/1.0e-4,1.0e-4/)
 real(r_kind) :: evpco=2.e-5 ! Zhao scheme evap coefficient for lg-scale rain
 real(r_kind) :: wminco(2)=(/1.e-5,1.e-5/) !  water and ice minimum threshold for Zhao 
 real(r_kind) :: flgmin(2)=(/0.2,0.2/) ! (Ferrier only) range of minimum large ice fraction
 logical :: old_monin = .false. ! flag for old Monin-Obhukov surface layer
 logical :: cnvgwd = .false. ! flag for convective gravity wave drag
 logical :: mom4ice = .false. ! flag for MOM4 sea-ice scheme
 logical :: shal_cnv = .true. ! use shallow convection?
 logical :: cal_pre = .false. ! true for huiya's precip type algorithm
 logical :: trans_trac = .true. ! convective transport of tracers? (RAS only)
 integer :: nst_fcst=0 ! 0 - AM only, 1 - uncoupled, 2 - coupled
 logical :: moist_adj = .false. 
 ! postphys=.false. means physics tendency computed at beginning of dynamics
 ! timestep, held constant in RK3 sub-steps.  If .true., tendency is applied
 ! after dynamics update as an adjustment (no physics tendencies in RK3
 ! sub-steps).
 ! postphys = .false. is similar to "process-split" physics, postphys=.true.
 ! is "time-split" physics (in the terminology of Williamson (2002): 
 ! http://journals.ametsoc.org/doi/abs/10.1175/1520-0493%282002%29130%3C2024%3ATSVPSC%3E2.0.CO%3B2).
 ! wrf uses postphys=.false. for everything except microphysics, which 
 ! is applied as an adjustment after the RK3 update (time-split).
 ! The operational GFS uses time split physics.
 ! time-split physics incurs the small extra cost of computing inverse transforms
 ! at the end of the dynamics time step.
 logical :: postphys = .true. 
 logical :: gloopb_filter = .true. ! apply spectral filter to physics tendencies

 real(r_kind) :: bkgd_vdif_m = 3.0 ! background vertical diffusion for momentum
 real(r_kind) :: bkgd_vdif_h = 1.0 ! background vertical diffusion for heat, q
 real(r_kind) :: bkgd_vdif_s = 0.2 ! sigma threshold for background mom. diffusn 
 ! if dt not given, but timestepsperhr is, dt=3600/timestepsperhr
 real(r_double) :: timestepsperhr = -1

 namelist/nam_mrf/initfile,sfcinitfile,fhmax,&
 deltim,dry,efold,ndiss,jablowill,heldsuarez,explicit,fhdfi,&
 fhout,fhzer,adiabatic,hdif_fac,hdif_fac2,fshk,ntrac,ntoz,ntclw,taustratdamp,&
 fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,polar_opt,&
 iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,psautco,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,mstrat,&
 pre_rad,bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,timestepsperhr,postphys,gloopb_filter

 contains

 subroutine read_namelist()
   integer lu,iret,ntracin
   logical idealized
   real(r_kind) tmax
   initfile=""
   sfcinitfile=""
   fhmax = 0
   fhout = 0
   fhzer = 0
   open(912,file='gfs_namelist',form='formatted')
   read(912,nam_mrf)
   close(912)
   if (initfile == "") then
      print *,'initfile must be specified in namelist'
      stop
   endif
   if (fhmax == 0) then
      print *,'fhmax must be specified in namelist'
      stop
   endif
   if (fhmax == 0) then
      print *,'fhout must be specified in namelist'
      stop
   endif
   if (deltim == 0) then
      if (timestepsperhr < 0) then
         print *,'deltim or timestepsperhr must be specified in namelist'
         stop
      else
         dt = 3600./timestepsperhr
      endif
   else
      dt=deltim
      timestepsperhr = nint(3600./dt)
   endif
   if (abs(3600.-dt*timestepsperhr) > 1.e-10) then
      print *,'1 hour must be an integer number of timesteps'
      stop
   endif
   if (postphys) then
      print *,'using time-split physics..'
   else
      print *,'using process-split physics..'
   endif
   lu = 7
   call sigio_sropen(lu,trim(initfile),iret)
   if (iret .ne. 0) then
      print *,'error opening ',trim(initfile),iret
      stop
   endif
   call sigio_srhead(lu,sighead,iret)
   if (iret .ne. 0) then
      print *,'error reading header from ',trim(initfile),iret
      stop
   else
      nlons = sighead%lonb
      nlats = sighead%latb
      nlevs = sighead%levs
      ntrunc = sighead%jcap
      tstart = sighead%fhour*3600.
      ntracin = sighead%ntrac
      idate_start = sighead%idate
      print *,'nlons,nlats,nlevs,ntrunc=',nlons,nlats,nlevs,ntrunc
      print *,'tstart=',tstart,' secs'
      print *,'idate_start=',idate_start
   endif 
   tmax = fhmax*3600. 
   ntmax = nint((tmax-tstart)/dt)
   ntout = fhout*timestepsperhr
   ntdfi = fhdfi*timestepsperhr
   print *,'output every ',ntout,' time steps'
   if (ntdfi > 0) print *,'digital filter half-window length',fhdfi,' hrs'
   idealized = jablowill .or. heldsuarez
   if (jablowill .and. heldsuarez) then
      print *,'conflicting namelist options'
      print *,'heldsuarez and jablowill both cannot be .true.'
      stop
   endif
   if (.not. idealized .and. (mod(ntdfi,int(fhswr*timestepsperhr)) .ne. 0 .or. &
       mod(ntdfi,int(fhlwr*timestepsperhr)) .ne. 0)) then
      print *,'middle of dfi window must be on a radiation time step'
      stop
   endif
   ! for these idealized tests, model is dry.
   if (jablowill) adiabatic = .true.
   if (jablowill) print *,'running jablonowsky and williamson test..'
   if (heldsuarez) print *,'running held-suarez test..'
   call sigio_sclose(lu,iret)
   ndimspec = (ntrunc+1)*(ntrunc+2)/2
   if (idealized .or. dry) ntrac=0
   if (.not. idealized .and. ntrac .ne. ntracin) then
     print *,ntracin,' tracers in input file, expecting',ntrac
     stop
   endif
   print *,'namelist nam_mrf:'
   write(6, nam_mrf)
 end subroutine read_namelist

end module params
