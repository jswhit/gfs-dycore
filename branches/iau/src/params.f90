module params
! holds model parameters
! Public subroutines:
! read_namelist: read namelist
 use kinds, only: r_kind,r_double
 use sigio_module, only: sigio_head, sigio_srhead, sigio_sropen, sigio_sclose

 implicit none
 private

 public :: read_namelist,initfile,sfcinitfile,fhmax,dt,ntmax,ndimspec,nlons,nlats,&
 tstart,ndiss,efold,nlevs,ntrunc,sighead,dry,explicit,heldsuarez,dcmip,&
 ntout,fhdfi,fhout,fhzer,idate_start,adiabatic,hdif_fac,hdif_fac2,fshk,ntrac,ntoz,ntclw,&
 pdryini,massfix,timestepsperhr,ncw,taustratdamp,polar_opt,ntdfi,gfsio_out,sigio_out,&
! gfs phys parameters.
 nmtvr,fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,&
 iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,psautco,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,mstrat,&
 pre_rad,bkgd_vdif_m,bkgd_vdif_s,bkgd_vdif_h,gloopb_filter,&
! iau parameters
 iau,iaufiles_fg,iaufiles_anl,iaufhrs,iau_delthrs,&
! vorticity confinement parameters
 vcamp,svc,svc_tau,svc_lscale,iseed_svc,svc_logit,&
! stochastic physics tendency parameters
 sppt,sppt_logit,sppt_tau,sppt_lscale,iseed_sppt, &
! additive stochastic humidity perturbations
 shum,shum_tau,shum_lscale,iseed_shum,clipsupersat
 character(len=500) :: initfile ! init cond filename
 character(len=500) :: sfcinitfile ! surface init cond filename
 integer            :: fhmax ! hours to run
 integer            :: fhout ! interval for IO
 integer            :: fhzer ! interval to zero accumulated arrays
! half window length (hrs) for digital filter launch (=0 mean no dfi)
 integer            :: fhdfi=0
 real(r_double)     :: deltim=0    ! namelist input time step (secs)
 real(r_double)     :: dt    ! time step (secs) (=deltim or 3600/timestepsperhr)
 real(r_kind) :: pdryini ! initial dry ps
 logical    :: massfix=.true. ! apply dry mass 'fixer'
 logical    :: gfsio_out=.false. ! write out 'gfsio' grib 1 files.
 logical    :: sigio_out=.true. ! write out 'sigma' spectral binary files.
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
 logical    :: iau = .false. ! iau forcing included
 integer    :: iau_delthrs = 6 ! iau time interval (to scale increments)
 character(len=120), dimension(7) ::  iaufiles_fg,iaufiles_anl
 real(r_kind), dimension(7) :: iaufhrs
 ! held-suarez forcing
 logical    :: heldsuarez = .false.
 ! dcmip test cases
 integer    :: dcmip = -1 ! 4x for baroclinic wave, 5x for tropical cyclone.
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
 real(r_kind) :: vcamp=0. ! vorticity confinement amplitude
 real(r_kind) :: svc=0.   ! stochastic vorticity confinement amplitude
 real(r_kind) :: svc_tau=0.      ! stochastic vorticity confinement time scale
 real(r_kind) :: svc_lscale=0.   ! stochastic vorticity confinement length scale
 integer :: iseed_svc=0 ! random seed for stochastic vc (zero means use clock)
 real(r_kind) :: sppt=0.  ! stochastic physics tendency amplitude
 real(r_kind) :: sppt_tau=0.  ! stochastic physics tendency time scale
 real(r_kind) :: sppt_lscale=0.  ! stochastic dynamics tendency length scale
 logical :: sppt_logit=.false. ! logit transform for sppt to bounded interval [-1,+1]
 logical :: svc_logit=.false.  ! logit transform for svc to bounded interval [-1,+1]
 integer :: iseed_sppt=0 ! random seed for sppt (0 means use system clock)
 integer :: iseed_shum=0 ! random seed for stochastic humid pert (0 means use system clock)
 real(r_kind) :: shum=0.  ! stochastic humidity pert amplitude
 real(r_kind) :: shum_tau=0.  ! stochastic humidity pert time scale
 real(r_kind) :: shum_lscale=0.  ! stochastic humidity pert length scale
 logical :: old_monin = .false. ! flag for old Monin-Obhukov surface layer
 logical :: cnvgwd = .false. ! flag for convective gravity wave drag
 logical :: mom4ice = .false. ! flag for MOM4 sea-ice scheme
 logical :: shal_cnv = .true. ! use shallow convection?
 logical :: cal_pre = .false. ! true for huiya's precip type algorithm
 logical :: trans_trac = .true. ! convective transport of tracers? (RAS only)
 integer :: nst_fcst=0 ! 0 - AM only, 1 - uncoupled, 2 - coupled
 logical :: moist_adj = .false. 
 logical :: gloopb_filter = .true. ! apply spectral filter to physics tendencies
 ! make sure stochastic perts don't create neg or supersat humidities.
 logical :: clipsupersat=.false. 

 real(r_kind) :: bkgd_vdif_m = 3.0 ! background vertical diffusion for momentum
 real(r_kind) :: bkgd_vdif_h = 1.0 ! background vertical diffusion for heat, q
 real(r_kind) :: bkgd_vdif_s = 0.2 ! sigma threshold for background mom. diffusn 
 ! if dt not given, but timestepsperhr is, dt=3600/timestepsperhr
 real(r_double) :: timestepsperhr = -1

 namelist/nam_mrf/initfile,sfcinitfile,fhmax,&
 massfix,deltim,dry,efold,ndiss,dcmip,heldsuarez,explicit,fhdfi,&
 fhout,fhzer,adiabatic,hdif_fac,hdif_fac2,fshk,ntrac,ntoz,ntclw,taustratdamp,&
 fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,polar_opt,&
 iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,psautco,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,mstrat,&
 pre_rad,bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,timestepsperhr,gloopb_filter,&
 vcamp,svc,svc_tau,svc_lscale,iseed_svc,sppt_tau,sppt,sppt_lscale,iseed_sppt,&
 clipsupersat,svc_logit,sppt_logit,shum,shum_tau,shum_lscale,iseed_shum,&
 gfsio_out,sigio_out,iau,iaufiles_fg,iaufiles_anl,iaufhrs,iau_delthrs

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
   iaufhrs = -1
   iaufiles_fg = ''
   iaufiles_anl = ''
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
   print *,'time step = ',dt,' seconds'
   if (abs(3600.-dt*timestepsperhr) > 1.e-10) then
      print *,'1 hour must be an integer number of timesteps'
      stop
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
      ! dry ps for mass fixer (if zero, compute in dyn_run)
      pdryini = sighead%pdryini*1000. ! convert to Pa from cb
      ! if pdryini is unrealistic, don't use it.
      if (pdryini .gt. 1.e5) then
          print *,'unrealistic pdryini in file ',pdryini,' resetting..'
          pdryini = 0.
       end if
      print *,'nlons,nlats,nlevs,ntrunc=',nlons,nlats,nlevs,ntrunc
      print *,'tstart=',tstart,' secs'
      print *,'pdryini=',pdryini
      print *,'idate_start=',idate_start
   endif 
   tmax = fhmax*3600. 
   ntmax = nint((tmax-tstart)/dt)
   ntout = fhout*timestepsperhr
   ntdfi = fhdfi*timestepsperhr
   print *,'output every ',ntout,' time steps'
   if (ntdfi > 0) print *,'digital filter half-window length',fhdfi,' hrs'
   idealized = dcmip >= 0 .or. heldsuarez
   if (dcmip >= 0 .and. heldsuarez) then
      print *,'conflicting namelist options'
      print *,'heldsuarez and dcmip both cannot be .true.'
      stop
   endif
   if (.not. idealized .and. (mod(ntdfi,int(fhswr*timestepsperhr)) .ne. 0 .or. &
       mod(ntdfi,int(fhlwr*timestepsperhr)) .ne. 0)) then
      print *,'middle of dfi window must be on a radiation time step'
      stop
   endif
   if (dcmip >= 0) print *,'running dcmip test...'
   if (heldsuarez) then
     print *,'running held-suarez test..'
     dry = .true.
   endif
   if (dcmip .eq. 41) dry = .true. ! dry baroclinic wave
   if (dcmip .eq. 42) dry = .false. ! dry baroclinic wave
   if (dcmip/10 .eq. 5) dry = .false. ! tropical cyclone
   call sigio_sclose(lu,iret)
   ndimspec = (ntrunc+1)*(ntrunc+2)/2
   if (dry) ntrac=0
   if (dcmip/10 .eq. 4) ntrac=3 ! passive tracers.
   if (dcmip    .eq. 51) ntrac=1 ! no passive tracers.
   if (.not. idealized .and. ntrac .ne. ntracin) then
     print *,ntracin,' tracers in input file, expecting',ntrac
     stop
   endif
   if (iau) then
     print *,'IAU forcing on'
   endif
   print *,'namelist nam_mrf:'
   write(6, nam_mrf)
 end subroutine read_namelist

end module params
