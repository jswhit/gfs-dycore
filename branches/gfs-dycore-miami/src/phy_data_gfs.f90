module phy_data
! data for GFS physics
! public subroutines:
! init_phydata: allocate and populate arrays.
! destroy_phydata: deallocate arrays.
! flx_init: initialize flx arrays.
 use kinds, only: r_kind, r_single, r_double
 use params, only: nlons,nlats,nlevs,ndimspec,sfcinitfile,nmtvr,ntoz,ntclw,num_p3d,num_p2d,&
      ntrunc,sighead,fhswr,fhlwr,idate_start,fhzer,dt,gloopb_filter
 use sfcio_module, only: sfcio_srohdc, sfcio_head, sfcio_data, sfcio_axdata, &
      sfcio_alhead, sfcio_aldata, sfcio_swohdc
 use physcons, only : tgice => con_tice, con_pi
 use shtns, only : lats,degree
 implicit none
 private
 public :: init_phydata, destroy_phydata, flx_init, wrtout_sfc, wrtout_flx

 type(sfcio_head),save :: sfchead
 integer, public :: lsoil, timeoz, pl_coeff, levozp, latsozp, &
                    levozc, latsozc, timeozc
 INTEGER, parameter, public :: MAX_SLOPETYP=30
 INTEGER, parameter, public :: MAX_SOILTYP=30
 INTEGER, parameter, public :: MAX_VEGTYP=30
 integer, parameter, public :: nfxr=33 ! dimension for rad fluxes array
 integer, parameter, public :: ncld=1  ! number of cloud types

 REAL, public :: SLOPE_DATA(MAX_SLOPETYP)
 REAL, public :: RSMTBL(MAX_VEGTYP)
 REAL, public :: RGLTBL(MAX_VEGTYP)
 REAL, public :: HSTBL(MAX_VEGTYP)
 REAL, public :: SNUPX(MAX_VEGTYP)
 REAL, public :: BB(MAX_SOILTYP)
 REAL, public :: DRYSMC(MAX_SOILTYP)
 REAL, public :: F11(MAX_SOILTYP)
 REAL, public :: MAXSMC(MAX_SOILTYP)
 REAL, public :: REFSMC(MAX_SOILTYP)
 REAL, public :: SATPSI(MAX_SOILTYP)
 REAL, public :: SATDK(MAX_SOILTYP)
 REAL, public :: SATDW(MAX_SOILTYP)
 REAL, public :: WLTSMC(MAX_SOILTYP)
 REAL, public :: QTZ(MAX_SOILTYP)
 LOGICAL, public :: LPARAM
 REAL, public :: ZBOT_DATA
 REAL, public :: SALP_DATA
 REAL, public :: CFACTR_DATA
 REAL, public :: CMCMAX_DATA
 REAL, public :: SBETA_DATA
 REAL, public :: RSMAX_DATA
 REAL, public :: TOPT_DATA
 REAL, public :: REFDK_DATA
 REAL, public :: FRZK_DATA
 INTEGER, public :: BARE
 INTEGER, public :: DEFINED_VEG
 INTEGER, public :: DEFINED_SOIL
 INTEGER, public :: DEFINED_SLOPE
 REAL, public :: FXEXP_DATA
 INTEGER, public :: NROOT_DATA(MAX_VEGTYP)
 REAL, public :: REFKDT_DATA
 REAL, public :: Z0_DATA(MAX_VEGTYP)
 REAL, public :: CZIL_DATA
 REAL, public :: LAI_DATA(MAX_VEGTYP)
 REAL, public :: CSOIL_DATA
 real(r_kind),public :: blatc, dphiozc ! ozone parameters (not used for prognostic ozone)

 ! ozone and mtnvar fields.
 real(r_kind),allocatable,public, dimension(:,:,:) :: hprime !mtnvar for grav wave drag
 real(r_kind),allocatable,public, dimension(:,:,:,:) :: ozplin !ozone P-L coeffs 
 real(r_kind),allocatable,public, dimension(:,:,:,:) :: phy_f3d !micro-physics 3d parameters
 real(r_kind),allocatable,public, dimension(:,:,:)   :: phy_f2d !micro-physics 2d parameters
 real(r_kind),allocatable,public, dimension(:) :: pl_lat,pl_pres,pl_time,ozddy !ozone lats,press,time
 real(r_kind),allocatable,public, dimension(:,:,:) :: cldcov ! 3d cloud fraction
 integer,allocatable,public,dimension(:) :: ozjindx1,ozjindx2 ! indices for lats
 !    solcon        - sun-earth distance adjusted solar constant (w/m2)  
 !    slag          - equation of time in radians                        
 !    sdec, cdec    - sin and cos of the solar declination angle         
 real(r_kind),public :: solcon,slag,sdec,cdec

 ! spectral filter for physics tendencies ("gloopb filter")
 real(r_kind), allocatable, public, dimension(:) :: bfilt
 ! 'sfc' fields (including orog).
 ! sea-surface temp K
 real(r_kind),allocatable,public, dimension(:,:)   :: tsea 
 ! soil volumetric water content (fraction)
 real(r_kind),allocatable,public, dimension(:,:,:) :: smc 
 ! snow depth in m
 real(r_kind),allocatable,public, dimension(:,:)   :: sheleg
 ! fractional snow cover 
 real(r_kind),allocatable,public, dimension(:,:)   :: sncovr
 ! soil temperature in K
 real(r_kind),allocatable,public, dimension(:,:,:) :: stc
 ! deep soil temperature in K
 real(r_kind),allocatable,public, dimension(:,:)   :: tg3
 ! roughness length in cm
 real(r_kind),allocatable,public, dimension(:,:)   :: zorl
 ! convective cloud cover (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: cv
 ! convective cloud bottom in kpa
 real(r_kind),allocatable,public, dimension(:,:)   :: cvb
 ! convective cloud top in kpa
 real(r_kind),allocatable,public, dimension(:,:)   :: cvt
 ! albedo for visible scattered (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: alvsf
 ! albedo for visible beam (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: alvwf
 ! albedo for near-IR scattered (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: alnsf
 ! albedo for near-IR beam (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: alnwf
 ! sea-land-ice mask (0-sea, 1-land, 2-ice)
 real(r_kind),allocatable,public, dimension(:,:)   :: slmsk
 ! vegetation fraction (fraction)
 real(r_kind),allocatable,public, dimension(:,:)   :: vfrac
 ! canopy water in m
 real(r_kind),allocatable,public, dimension(:,:)   :: canopy
 ! 10-meter wind speed over lowest model wind speed
 real(r_kind),allocatable,public, dimension(:,:)   :: f10m
 ! 2-meter temperature in K
 real(r_kind),allocatable,public, dimension(:,:)   :: t2m
 ! 2-meter specific humidity in kg/kg
 real(r_kind),allocatable,public, dimension(:,:)   :: q2m
 ! vegetation type in integer 1-13
 real(r_kind),allocatable,public, dimension(:,:)   :: vtype
 ! soil type in integer 1-9
 real(r_kind),allocatable,public, dimension(:,:)   :: stype
 ! fractional coverage with strong cosz dependency
 ! (multiplies alvsf and alnsf in albedo calculation in
 ! setalb/radiation_surface.f90)
 real(r_kind),allocatable,public, dimension(:,:)   :: facsf
 ! fractional coverage with weak cosz dependency  
 ! (multiplies alvwf and alnwf in albedo calculation in 
 ! setalb/radiation_surface.f90)
 real(r_kind),allocatable,public, dimension(:,:)   :: facwf
 ! surface layer friction velocity
 real(r_kind),allocatable,public, dimension(:,:)   :: uustar
 ! F_m parameter from PBL scheme = LOG((Z0MAX(I)+Z1(I)) / Z0MAX(I))
 real(r_kind),allocatable,public, dimension(:,:)   :: ffmm
 ! F_h parameter from PBL scheme = LOG((ZTMAX(I)+Z1(I)) / ZTMAX(I))
 real(r_kind),allocatable,public, dimension(:,:)   :: ffhh
 ! sea-ice thickness
 real(r_kind),allocatable,public, dimension(:,:)   :: hice
 ! sea-ice concentration
 real(r_kind),allocatable,public, dimension(:,:)   :: fice
 ! sea-ice temperature
 real(r_kind),allocatable,public, dimension(:,:)   :: tisfc
 ! accumulated total precipitation (kg/m2) 
 real(r_kind),allocatable,public, dimension(:,:)   :: tprcp
 ! snow/rain flag for precipitation (0 rain, 1 snow)
 real(r_kind),allocatable,public, dimension(:,:)   :: srflag
 ! actual snow depth (mm) over land/sea ice 
 real(r_kind),allocatable,public, dimension(:,:)   :: snwdph
 ! liquid soil moisture content (fraction)
 real(r_kind),allocatable,public, dimension(:,:,:) :: slc
 ! minimum areal coverage of green veg
 real(r_kind),allocatable,public, dimension(:,:)   :: shdmin
 ! maximum areal coverage of green veg
 real(r_kind),allocatable,public, dimension(:,:)   :: shdmax
 ! integer class of sfc slope 
 real(r_kind),allocatable,public, dimension(:,:)   :: slope
 ! maximum (deep) snow albedo 
 real(r_kind),allocatable,public, dimension(:,:)   :: snoalb
 ! spectrally filtered orography in m
 real(r_kind),allocatable,public, dimension(:,:)   :: oro
 ! unfiltered orography in m
 real(r_kind),allocatable,public, dimension(:,:)   :: oro_uf

 ! 'flx' fields.
 ! snow-free albedo
 real(r_kind),public,allocatable:: SFALB(:,:)
 ! total sky surface downward sw flux in w/m**2
 real(r_kind),public,allocatable:: SFCDSW(:,:)
 ! mean cos of zenith angle over rad call period 
 real(r_kind),public,allocatable:: COSZEN(:,:)
 ! precipitable water
 real(r_kind),public,allocatable:: PWAT(:,:)
 ! minimum 2m temp
 real(r_kind),public,allocatable:: TMPMIN(:,:)
 ! maximum 2m temp
 real(r_kind),public,allocatable:: TMPMAX(:,:)
 ! minimum 2m specific humidity
 real(r_kind),public,allocatable:: SPFHMIN(:,:)
 ! maximum 2m specific humidity
 real(r_kind),public,allocatable:: SPFHMAX(:,:)
 ! u component of surface stress
 real(r_kind),public,allocatable:: DUSFC(:,:)
 ! v component of surface stress
 real(r_kind),public,allocatable:: DVSFC(:,:)
 ! sensible heat flux (w/m2) 
 real(r_kind),public,allocatable:: DTSFC(:,:)
 ! latent heat flux (w/m2)  
 real(r_kind),public,allocatable:: DQSFC(:,:)
 ! time accumulated sfc dn lw flux ( w/m**2 )
 real(r_kind),public,allocatable:: DLWSFC(:,:)
 ! time accumulated sfc up lw flux ( w/m**2 ) 
 real(r_kind),public,allocatable:: ULWSFC(:,:)
 ! ground conductive heat flux 
 real(r_kind),public,allocatable:: GFLUX(:,:)
 ! total water runoff 
 real(r_kind),public,allocatable:: RUNOFF(:,:)
 ! potential evaporation
 real(r_kind),public,allocatable:: EP(:,:)
 ! cloud workfunction (valid only with sas)
 real(r_kind),public,allocatable:: CLDWRK(:,:)
 ! vertically integrated u change by OGWD  
 real(r_kind),public,allocatable:: DUGWD(:,:)
 ! vertically integrated v change by OGWD
 real(r_kind),public,allocatable:: DVGWD(:,:)
 ! surface pressure (kPa) average over physics time step
 real(r_kind),public,allocatable:: PSMEAN(:,:)
 ! total precipitation rate
 real(r_kind),public,allocatable:: GESHEM(:,:)
 ! convective precipitation rate
 real(r_kind),public,allocatable:: BENGSH(:,:)
 ! total sky sfc net sw flux into ground in w/m**2
 real(r_kind),public,allocatable:: SFCNSW(:,:)
 ! total sky surface downward lw flux in w/m**2 
 real(r_kind),public,allocatable:: SFCDLW(:,:)
 ! surface air temp during lw calculation in k
 real(r_kind),public,allocatable:: TSFLW(:,:)
 ! instantaneous surface pressure (kPa)
 real(r_kind),public,allocatable:: PSURF(:,:)
 ! 10m u
 real(r_kind),public,allocatable:: U10M(:,:)
 ! 10m v
 real(r_kind),public,allocatable:: V10M(:,:)
 ! pbl height
 real(r_kind),public,allocatable:: HPBL(:,:)
 ! thermal exchange coefficient
 real(r_kind),public,allocatable:: CHH(:,:)
 !  momentum exchange coefficient
 real(r_kind),public,allocatable:: CMM(:,:)
 ! instantaneous sfc potential evaporation
 real(r_kind),public,allocatable:: EPI(:,:)
 ! instantaneous sfc dnwd lw flux ( w/m**2 )
 real(r_kind),public,allocatable:: DLWSFCI(:,:)
 ! instantaneous sfc upwd lw flux ( w/m**2 )
 real(r_kind),public,allocatable:: ULWSFCI(:,:)
 ! instantaneous sfc upwd sw flux ( w/m**2 )
 real(r_kind),public,allocatable:: USWSFCI(:,:)
 ! instantaneous sfc dnwd sw flux ( w/m**2 ) 
 real(r_kind),public,allocatable:: DSWSFCI(:,:)
 ! instantaneous sfc sensible heat flux 
 real(r_kind),public,allocatable:: DTSFCI(:,:)
 ! instantaneous sfc latent heat flux 
 real(r_kind),public,allocatable:: DQSFCI(:,:)
 ! instantaneous sfc ground heat flux
 real(r_kind),public,allocatable:: GFLUXI(:,:)
 ! surface water runoff (from lsm) 
 real(r_kind),public,allocatable:: SRUNOFF(:,:)
 ! layer 1 temperature (K) 
 real(r_kind),public,allocatable:: T1(:,:)
 ! layer 1 specific humidity (kg/kg)
 real(r_kind),public,allocatable:: Q1(:,:)
 ! layer 1 zonal wind (m/s) 
 real(r_kind),public,allocatable:: U1(:,:)
 ! layer 1 merdional wind (m/s)
 real(r_kind),public,allocatable:: V1(:,:)
 ! layer 1 height (m) 
 real(r_kind),public,allocatable:: ZLVL(:,:)
 ! Direct evaporation from bare soil
 real(r_kind),public,allocatable:: EVBSA(:,:)
 ! Canopy water evaporation
 real(r_kind),public,allocatable:: EVCWA(:,:)
 ! Transpiration land surface 
 real(r_kind),public,allocatable:: TRANSA(:,:)
 ! Snow Sublimation land surface  
 real(r_kind),public,allocatable:: SBSNOA(:,:)
 ! Snow Cover (fraction) 
 real(r_kind),public,allocatable:: SNOWCA(:,:)
 ! Snow phase-change heat fluxi  
 real(r_kind),public,allocatable:: SNOHFA(:,:)
 ! total column soil moisture
 real(r_kind),public,allocatable:: SOILM(:,:)
 ! wilting point (volumetric) 
 real(r_kind),public,allocatable:: SMCWLT2(:,:)
 ! soil moisture threshold (volumetric)
 real(r_kind),public,allocatable:: SMCREF2(:,:)
 ! sunshine duration time
 real(r_kind),public,allocatable:: suntim(:,:)       
 ! surface emissivity
 real(r_kind),public,allocatable:: sfcemis(:,:)      
 ! Average VOL soil moist content(frac) layer 10cm -> 0cm
 real(r_kind),public,allocatable:: gsoil(:,:)
 ! Average temperature at 2 meter (K) 
 real(r_kind),public,allocatable:: gtmp2m(:,:)
 ! Average Frictional Velocity (m/s)
 real(r_kind),public,allocatable:: gustar(:,:)
 ! Average PBL height
 real(r_kind),public,allocatable:: gpblh(:,:)
 ! Average 10m u
 real(r_kind),public,allocatable:: gu10m(:,:)
 ! Average 10m v
 real(r_kind),public,allocatable:: gv10m(:,:)
 ! Average Surface roughness (m)
 real(r_kind),public,allocatable:: gzorl(:,:)
 ! Average Land Sea surface (fraction)
 real(r_kind),public,allocatable:: goro(:,:)
 ! radiative fluxes and heating (need to be saved between radiation time steps)
 real(r_kind),public,allocatable:: fluxr(:,:,:)
 real(r_kind),public,allocatable:: swh(:,:,:)
 real(r_kind),public,allocatable:: hlw(:,:,:)
 logical :: is_allocated = .false.

 contains

 subroutine init_phydata()
   type(sfcio_data) data
   integer, parameter ::  nread=14
   integer, parameter ::  nmtn=24
   integer, parameter ::  kozpl=28
   integer iret,vegtyp,nf0,nf1
   real(r_kind) rsnow,lat,fd2
   real(r_single), allocatable, dimension(:,:,:) :: hprime4
   real(r_single), allocatable, dimension(:) :: tempin,pl_lat4,pl_pres4,pl_time4
   integer i,j,k,n,nc

   if (sfcinitfile == "") then
      print *,'sfcinitfile must be specified in namelist'
      stop
   endif
   print *,' nread=',nread,' sfcinitfile=',trim(sfcinitfile)
   call sfcio_srohdc(nread,sfcinitfile,sfchead,data,iret)
   write(6,99) nread,sfchead%fhour,sfchead%idate,&
   sfchead%lonb,sfchead%latb,sfchead%lsoil,sfchead%ivs,iret
99 FORMAT(1H ,'in fixio nread=',i3,2x,'HOUR=',f8.2,3x,'IDATE=', &
   4(1X,I4),4x,'lonsfc,latsfc,lsoil,ivssfc,iret=',5i8)
   if(iret.ne.0) then
     write(6,*) 'error reading ',trim(sfcinitfile)
     stop
   endif
   if (sfchead%lonb .ne. nlons .or. sfchead%latb .ne. nlats) then
     write(6,*) 'sfcfile dims',sfchead%lonb,sfchead%latb,' expecting ',nlons,nlats
     stop
   endif
   lsoil = sfchead%lsoil

   ! only do this part on first call when arrays not allocated.
   ! when called after dfi, it will skip this bloack
   if (.not. is_allocated) then

   is_allocated = .true.
   ! spectral filter for physics tendencies ("gloopb filter")
   allocate(bfilt(ndimspec))
   bfilt=1 ! (no filtering)
   if (gloopb_filter) then
      nf0 = (ntrunc+1)*2/3  ! highest wavenumber gloopb filter keeps fully
      !nf0 = (ntrunc+1)*9/10 ! used for semi-lagrangian
      nf1 = (ntrunc+1)      ! lowest wavenumber gloopb filter removes fully
      fd2 = 1./(nf1-nf0)**2
      do nc=1,ndimspec
         n = degree(nc)
         bfilt(nc) = max(1.-fd2*max(n-nf0,0)**2,0.)
      enddo
   endif

   ! set land model parameters.
   call set_soilveg()

   ! microphysics arrays
   allocate(phy_f3d(nlons,nlats,nlevs,num_p3d))
   allocate(phy_f2d(nlons,nlats,num_p2d))

   ! 3d cloud fraction diagnosed by radiation.
   allocate(cldcov(nlons,nlats,nlevs))

   ! read mtnvar
   allocate(hprime(nlons,nlats,nmtvr))
   allocate(hprime4(nlons,nlats,nmtvr))
   open(nmtn,file='mtnvar.dat',form='unformatted',convert='big_endian')
   READ(nmtn) hprime4
   close(nmtn)
   hprime = hprime4
   deallocate(hprime4)

   ! read ozone forcing for prognostic ozone
   if (ntoz > 0) then ! prognostic ozone
      open(kozpl,file='o3forcing.dat',form='unformatted',convert='big_endian')
      rewind (kozpl)
      read (kozpl) pl_coeff, latsozp, levozp, timeoz
      allocate (pl_lat(latsozp), pl_pres(levozp),pl_time(timeoz+1))
      allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1))
      rewind (kozpl)
      read (kozpl) pl_coeff, latsozp, levozp, timeoz, pl_lat4, pl_pres4,  &
                   pl_time4
      pl_pres(:) = log(100.0*pl_pres4(:))  ! Natural log of pres in Pa
      pl_lat(:)  = pl_lat4(:)
      pl_time(:) = pl_time4(:)
      ! these are not relevant for prognostic ozone, but some 
      ! physics modules import them.
      blatc   = 0.0; dphiozc = 0.0
      levozc  = 17
      latsozc = 18
      timeozc = 1
      allocate(tempin(latsozp))
      allocate(ozplin(latsozp,levozp,pl_coeff,timeoz))
      DO I=1,timeoz
        do n=1,pl_coeff
          DO k=1,levozp
            READ(kozpl) tempin
            ozplin(:,k,n,i) = tempin(:)
          ENDDO
        enddo
      ENDDO
      deallocate(tempin)
      allocate(ozjindx1(nlats),ozjindx2(nlats),ozddy(nlats))
      do j=1,nlats
         lat = (180./con_pi)*lats(1,j)
         ozjindx2(j) = latsozp + 1
         do i=1,latsozp
            if (lat .lt. pl_lat(i)) then
              ozjindx2(j) = i
              exit
            endif
         enddo
         ozjindx1(j) = max(ozjindx2(j)-1,1)
         ozjindx2(j) = min(ozjindx2(j),latsozp)
         if (ozjindx2(j) .ne. ozjindx1(j)) then
           OZDDY(j) = (lat             - pl_lat(ozjindx1(j))) &
                  / (pl_lat(ozjindx2(j)) - pl_lat(ozjindx1(j)))
         else
           ozddy(j) = 1.0
         endif
      enddo
      close(kozpl)
    end if
!   reading the grib orography
    allocate(oro(nlons,nlats))
    allocate(oro_uf(nlons,nlats))
    CALL ORORD(101,nlons,nlats,oro,'orography')
    print *,'read grb orography'
!   read unfiltered orography
    CALL ORORD(101,nlons,nlats,oro_uf,'orography_uf')
    print *,'read grb orography_uf'
    ! allocate arrays for data in sfc file.
    allocate(tsea(nlons,nlats))
    allocate(smc(nlons,nlats,lsoil))
    allocate(stc(nlons,nlats,lsoil))
    allocate(tg3(nlons,nlats))
    allocate(zorl(nlons,nlats))
    allocate(cv(nlons,nlats))
    allocate(cvb(nlons,nlats))
    allocate(cvt(nlons,nlats))
    allocate(alvsf(nlons,nlats))
    allocate(alvwf(nlons,nlats))
    allocate(alnwf(nlons,nlats))
    allocate(alnsf(nlons,nlats))
    allocate(slmsk(nlons,nlats))
    allocate(vfrac(nlons,nlats))
    allocate(canopy(nlons,nlats))
    allocate(f10m(nlons,nlats))
    allocate(t2m(nlons,nlats))
    allocate(q2m(nlons,nlats))
    allocate(vtype(nlons,nlats))
    allocate(stype(nlons,nlats))
    allocate(facsf(nlons,nlats))
    allocate(facwf(nlons,nlats))
    allocate(uustar(nlons,nlats))
    allocate(ffmm(nlons,nlats))
    allocate(ffhh(nlons,nlats))
    allocate(hice(nlons,nlats))
    allocate(fice(nlons,nlats))
    allocate(tisfc(nlons,nlats))
    allocate(tprcp(nlons,nlats))
    allocate(srflag(nlons,nlats))
    allocate(snwdph(nlons,nlats))
    allocate(slc(nlons,nlats,lsoil))
    allocate(shdmin(nlons,nlats))
    allocate(shdmax(nlons,nlats))
    allocate(slope(nlons,nlats))
    allocate(snoalb(nlons,nlats))
    allocate(sheleg(nlons,nlats))
    allocate(sncovr(nlons,nlats))
    allocate(               &
    SFALB   (nlons,nlats), &
    SFCDSW  (nlons,nlats), &
    COSZEN  (nlons,nlats), &
    PWAT    (nlons,nlats), &
    TMPMIN  (nlons,nlats), &
    TMPMAX  (nlons,nlats), &
    SPFHMIN (nlons,nlats), &
    SPFHMAX (nlons,nlats), &
    DUSFC   (nlons,nlats), &
    DVSFC   (nlons,nlats), &
    DTSFC   (nlons,nlats), &
    DQSFC   (nlons,nlats), &
    DLWSFC  (nlons,nlats), &
    ULWSFC  (nlons,nlats), &
    GFLUX   (nlons,nlats), &
    RUNOFF  (nlons,nlats), &
    EP      (nlons,nlats), &
    CLDWRK  (nlons,nlats), &
    DUGWD   (nlons,nlats), &
    DVGWD   (nlons,nlats), &
    PSMEAN  (nlons,nlats), &
    GESHEM  (nlons,nlats), &
    BENGSH  (nlons,nlats), &
    SFCNSW  (nlons,nlats), &
    SFCDLW  (nlons,nlats), &
    TSFLW   (nlons,nlats), &
    PSURF   (nlons,nlats), &
    U10M    (nlons,nlats), &
    V10M    (nlons,nlats), &
    HPBL    (nlons,nlats), &
    CHH     (nlons,nlats), &
    CMM     (nlons,nlats), &
    EPI     (nlons,nlats), &
    DLWSFCI (nlons,nlats), &
    ULWSFCI (nlons,nlats), &
    USWSFCI (nlons,nlats), &
    DSWSFCI (nlons,nlats), &
    DTSFCI  (nlons,nlats), &
    DQSFCI  (nlons,nlats), &
    GFLUXI  (nlons,nlats), &
    SRUNOFF (nlons,nlats), &
    T1      (nlons,nlats), &
    Q1      (nlons,nlats), &
    U1      (nlons,nlats), &
    V1      (nlons,nlats), &
    ZLVL    (nlons,nlats), &
    EVBSA   (nlons,nlats), &
    EVCWA   (nlons,nlats), &
    TRANSA  (nlons,nlats), &
    SBSNOA  (nlons,nlats), &
    SNOWCA  (nlons,nlats), &
    SOILM   (nlons,nlats), &
    SNOHFA  (nlons,nlats), &
    SMCWLT2 (nlons,nlats), &
    SMCREF2 (nlons,nlats), &
    suntim  (nlons,nlats), &                
    sfcemis (nlons,nlats), &                
    gsoil   (nlons,nlats), &
    gtmp2m  (nlons,nlats), &
    gustar  (nlons,nlats), &
    gpblh   (nlons,nlats), &
    gu10m   (nlons,nlats), &
    gv10m   (nlons,nlats), &
    gzorl   (nlons,nlats), &
    goro    (nlons,nlats), &
    stat=iret)
    allocate(fluxr(nlons,nlats,nfxr))
    allocate(swh(nlons,nlats,nlevs))
    allocate(hlw(nlons,nlats,nlevs))

   endif ! skip to here if arrays already allocated

   phy_f3d=0; phy_f2d=0
   cldcov=0.

   tsea = data%tsea
   smc = data%smc
   stc = data%stc
   tg3 = data%tg3
   zorl = data%zorl
   !cv = data%cv
   !cvb = data%cvb
   !cvt = data%cvt
   ! not needed anyway (only used for diagnostic clouds).
   ! file contains missing values
   cv = 0; cvb = 0; cvt = 0;
   alvsf = data%alvsf
   alvwf = data%alvwf
   alnwf = data%alnwf
   alnsf = data%alnsf
   slmsk = data%slmsk
   vfrac = data%vfrac
   canopy = data%canopy
   f10m = data%f10m
   t2m = data%t2m
   q2m = data%q2m
   vtype = data%vtype
   stype = data%stype
   facsf = data%facsf
   facwf = data%facwf
   uustar = data%uustar
   ffmm = data%ffmm
   ffhh = data%ffhh
   hice = data%hice
   fice = data%fice
   tisfc = data%tisfc

   ! if tisfc<0, determine from sst, ice concentration and tgice (constant)
   if (tisfc(1,1) < 0.) then
      do j=1,nlats
      do i=1,nlons
         tisfc(i,j) = tsea(i,j)
         if (slmsk(i,j) >= 2. .and. fice(i,j) >= 0.5) then
            tisfc(i,j) = tsea(i,j) - tgice*(1.-fice(i,j))/fice(i,j)
            tisfc(i,j) = min(tisfc(i,j),tgice)
         end if
      enddo
      enddo
   endif


   tprcp = data%tprcp
   srflag = data%srflag
   snwdph = data%snwdph
   slc = data%slc
   shdmin = data%shdmin
   shdmax = data%shdmax
   slope = data%slope
   snoalb = data%snoalb
   ! orog read directly from grib file.
   !oro = data%oro
   sheleg = data%sheleg

! initialize snow fraction(sheleg is in mm)
   do j=1,nlats
   do i=1,nlons
      sncovr(i,j) = 0.!
      if (slmsk(i,j) > 0.001) then
         vegtyp = VTYPE(i,j)
         if( vegtyp .eq. 0 ) vegtyp = 7	
         RSNOW  = 0.001*SHELEG(i,j)/SNUPX(vegtyp)
         IF (0.001*SHELEG(i,j) < SNUPX(vegtyp)) THEN
           SNCOVR(i,j) = 1.0 - ( EXP(-SALP_DATA*RSNOW) - &
                                 RSNOW*EXP(-SALP_DATA))
         ELSE
           SNCOVR(i,j) = 1.0
         ENDIF
      endif
   enddo
   enddo

   call sfcio_axdata(data,iret)

   call flx_init()

 end subroutine init_phydata

 subroutine flx_init()
    TMPMIN  = 1.e10
    TMPMAX  = 0.
    SPFHMIN = 1.e10
    SPFHMAX = 0.
    GESHEM  = 0.
    BENGSH  = 0.
    DUSFC   = 0.
    DVSFC   = 0.
    DTSFC   = 0.
    DQSFC   = 0.
    DLWSFC  = 0.
    ULWSFC  = 0.
    suntim  = 0.
    GFLUX   = 0.
    RUNOFF  = 0.
    EP      = 0.
    CLDWRK  = 0.
    DUGWD   = 0.
    DVGWD   = 0.
    PSMEAN  = 0.
    EVBSA   = 0.
    EVCWA   = 0.
    TRANSA  = 0.
    SBSNOA  = 0.
    SNOWCA  = 0.
    SRUNOFF = 0.
    SNOHFA  = 0.
    gsoil   = 0.
    gtmp2m  = 0.
    gustar  = 0.
    gpblh   = 0.
    gu10m   = 0.
    gv10m   = 0.
    gzorl   = 0.
    goro    = 0.
    fluxr   = 0.
    swh     = 0.
    hlw     = 0.
    cv      = 0.
    cvb     = 0.
    cvt     = 0.
 end subroutine flx_init

 subroutine set_soilveg()
   integer i
   REAL WLTSMC1,REFSMC1
! ----------------------------------------------------------------------
! SET TWO SOIL MOISTURE WILT, SOIL MOISTURE REFERENCE PARAMETERS
! ----------------------------------------------------------------------
   REAL SMLOW
   REAL SMLOW_DATA
   DATA SMLOW_DATA /0.5/

   REAL SMHIGH
   REAL SMHIGH_DATA
   DATA SMHIGH_DATA /6.0/
   NAMELIST /SOIL_VEG/ SLOPE_DATA, RSMTBL, RGLTBL, HSTBL, SNUPX, &
    BB, DRYSMC, F11, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, &
    WLTSMC, QTZ, LPARAM, ZBOT_DATA, SALP_DATA, CFACTR_DATA, &
    CMCMAX_DATA, SBETA_DATA, RSMAX_DATA, TOPT_DATA, &
    REFDK_DATA, FRZK_DATA, BARE, DEFINED_VEG, DEFINED_SOIL, &
    DEFINED_SLOPE, FXEXP_DATA, NROOT_DATA, REFKDT_DATA, Z0_DATA, &
    CZIL_DATA, LAI_DATA, CSOIL_DATA 

   SLOPE_DATA =(/0.1,  0.6, 1.0, 0.35, 0.55, 0.8, &
   	       0.63, 0.0, 0.0, 0.0,  0.0,  0.0, &
   	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0, &
   	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0, &
   	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0/)
!  RSMTBL =(/300.0, 175.0, 175.0, 300.0,  70.0, 70.0, &
!           20.0,  70.0,  70.0,  70.0,  70.0, 20.0, &
!           70.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
!            0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
!            0.0,   0.0,   0.0,   0.0,   0.0,  0.0/)
! Sep 5 2012 bugfix (svn r 20728)
   RSMTBL =(/300.0, 175.0, 175.0, 300.0, 300.0, 70.0, &
             45.0, 225.0, 225.0, 225.0, 400.0, 45.0, &
            150.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
              0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
              0.0,   0.0,   0.0,   0.0,   0.0,  0.0/)
   RGLTBL =(/30.0,  30.0,  30.0,  30.0,  30.0,  65.0, &
    	  100.0, 100.0, 100.0, 100.0, 100.0, 100.0, &
    	  100.0,   0.0,   0.0,   0.0,	0.0,   0.0, &
    	    0.0,   0.0,   0.0,   0.0,	0.0,   0.0, &
    	    0.0,   0.0,   0.0,   0.0,	0.0,   0.0/)
   HSTBL =(/41.69, 54.53, 51.93, 47.35,  47.35, 54.53, &
    	  36.35, 42.00, 42.00, 42.00,  42.00, 36.35, &
    	  42.00,  0.00,  0.00,  0.00,	0.00,  0.00, &
    	   0.00,  0.00,  0.00,  0.00,	0.00,  0.00, &
    	   0.00,  0.00,  0.00,  0.00,	0.00,  0.00/)
   SNUPX  =(/0.040, 0.040, 0.040, 0.040, 0.040, 0.040, &
        0.020, 0.020, 0.020, 0.020, 0.013, 0.020, &
        0.013, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	   0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	   0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
   BB	 =(/4.26,  8.72, 11.55, 4.74, 10.73,  8.17, &
    	  6.77,  5.25,  4.26, 0.00,  0.00,  0.00, &
    	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00, &
    	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00, &
    	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00/)
   DRYSMC=(/0.029, 0.119, 0.139, 0.047, 0.100, 0.103, &
    	  0.069, 0.066, 0.029, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
   F11  =(/-0.999, -1.116, -2.137, -0.572, -3.201, -1.302, &
    	 -1.519, -0.329, -0.999,  0.000,  0.000,  0.000, &
    	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, &
    	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, &
    	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000/)
   MAXSMC=(/0.421, 0.464, 0.468, 0.434, 0.406, 0.465, &
    	  0.404, 0.439, 0.421, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
!
! ----------------------------------------------------------------------
! THE FOLLOWING 5 PARAMETERS ARE DERIVED LATER IN REDPRM.F FROM THE SOIL
! DATA, AND ARE JUST GIVEN HERE FOR REFERENCE AND TO FORCE STATIC
! STORAGE ALLOCATION. -DAG LOHMANN, FEB. 2001
! ----------------------------------------------------------------------
   REFSMC=(/0.248, 0.368, 0.398, 0.281, 0.321, 0.361, &
       0.293, 0.301, 0.248, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
! ---------------------------------------------------------------------
! SOIL TEXTURE-RELATED ARRAYS.
! ----------------------------------------------------------------------
   SATPSI=(/0.04, 0.62, 0.47, 0.14, 0.10, 0.26, &
    	  0.14, 0.36, 0.04, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00/)
   SATDK =(/1.41E-5, 0.20E-5, 0.10E-5, 0.52E-5, 0.72E-5, &
    	  0.25E-5, 0.45E-5, 0.34E-5, 1.41E-5, 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00/)
   QTZ   =(/0.82, 0.10, 0.25, 0.60, 0.52, 0.35, &
    	  0.60, 0.40, 0.82, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
    	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00/)
   WLTSMC=(/0.029, 0.119, 0.139, 0.047, 0.100, 0.103, &
    	  0.069, 0.066, 0.029, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
   SATDW =(/5.71E-6, 2.33E-5, 1.16E-5, 7.95E-6, 1.90E-5, &
    	  1.14E-5, 1.06E-5, 1.46E-5, 5.71E-6, 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00, &
    	  0.00   , 0.00   , 0.00   , 0.00   , 0.00/)

   LPARAM =.TRUE.
   ZBOT_DATA =-8.0
   SALP_DATA =4.0
   CFACTR_DATA =0.5
   CMCMAX_DATA =0.5E-3
   SBETA_DATA =-2.0
   RSMAX_DATA =5000.0
   TOPT_DATA =298.0
   REFDK_DATA =2.0E-6
   FRZK_DATA =0.15
   BARE =11

! ----------------------------------------------------------------------
! NUMBER OF DEFINED SOIL-, VEG-, AND SLOPETYPS USED.
! ----------------------------------------------------------------------

   DEFINED_VEG=13
   DEFINED_SOIL=9
   DEFINED_SLOPE=9
 
   FXEXP_DATA =2.0
   !NROOT_DATA =(/4,4,4,4,4,4,4,4,4,4,4,4,4,0,0, &
   ! Sep 5 2012 bugfix (svn r 20728)
   NROOT_DATA =(/4,4,4,4,4,4,3,3,3,2,3,3,2,0,0, &
   	        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
   REFKDT_DATA =3.0
! ----------------------------------------------------------------------
! VEGETATION CLASS-RELATED ARRAYS
! ----------------------------------------------------------------------
   Z0_DATA =(/2.653, 0.826, 0.563, 1.089, 0.854, 0.856, &
    	    0.035, 0.238, 0.065, 0.076, 0.011, 0.035, &
    	    0.011, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
    	    0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
   CZIL_DATA =0.075
   LAI_DATA =(/3.0, 3.0, 3.0, 3.0, 3.0, 3.0, &
    	     3.0, 3.0, 3.0, 3.0, 3.0, 3.0, &
    	     3.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    	     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    	     0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
   CSOIL_DATA =2.00E+6
! ----------------------------------------------------------------------
! READ NAMELIST FILE TO OVERRIDE DEFAULT PARAMETERS ONLY ONCE.
! ----------------------------------------------------------------------
   WRITE(6,*) 'READ NAMELIST SOIL_VEG'
   open(912,form="formatted",file="gfs_namelist")
   READ(912, SOIL_VEG)
   close(912)
   IF (DEFINED_SOIL .GT. MAX_SOILTYP) THEN
      WRITE(6,*) 'Warning: DEFINED_SOIL too large in namelist'
      STOP 222
   ENDIF
   IF (DEFINED_VEG .GT. MAX_VEGTYP) THEN
      WRITE(6,*) 'Warning: DEFINED_VEG too large in namelist'
      STOP 222
   ENDIF
   IF (DEFINED_SLOPE .GT. MAX_SLOPETYP) THEN
      WRITE(6,*) 'Warning: DEFINED_SLOPE too large in namelist'
      STOP 222
   ENDIF
     
   SMLOW = SMLOW_DATA
   SMHIGH = SMHIGH_DATA
      
   DO I = 1,DEFINED_SOIL
      SATDW(I)  = BB(I)*SATDK(I)*(SATPSI(I)/MAXSMC(I))
      F11(I) = ALOG10(SATPSI(I)) + BB(I)*ALOG10(MAXSMC(I)) + 2.0
      REFSMC1 = MAXSMC(I)*(5.79E-9/SATDK(I)) &
             **(1.0/(2.0*BB(I)+3.0))
      REFSMC(I) = REFSMC1 + (MAXSMC(I)-REFSMC1) / SMHIGH
      WLTSMC1 = MAXSMC(I) * (200.0/SATPSI(I))**(-1.0/BB(I))
      WLTSMC(I) = WLTSMC1 - SMLOW * WLTSMC1
         
!  ----------------------------------------------------------------------
!  CURRENT VERSION DRYSMC VALUES THAT EQUATE TO WLTSMC.
!  FUTURE VERSION COULD LET DRYSMC BE INDEPENDENTLY SET VIA NAMELIST.
!  ----------------------------------------------------------------------
      DRYSMC(I) = WLTSMC(I)
   END DO
         
 end subroutine set_soilveg

 subroutine destroy_phydata()
   deallocate(tsea)
   deallocate(smc)
   deallocate(stc)
   deallocate(tg3)
   deallocate(zorl)
   deallocate(cv)
   deallocate(cvb)
   deallocate(cvt)
   deallocate(alvsf)
   deallocate(alvwf)
   deallocate(alnwf)
   deallocate(alnsf)
   deallocate(slmsk)
   deallocate(vfrac)
   deallocate(canopy)
   deallocate(f10m)
   deallocate(t2m)
   deallocate(q2m)
   deallocate(vtype)
   deallocate(stype)
   deallocate(facsf)
   deallocate(facwf)
   deallocate(uustar)
   deallocate(ffmm)
   deallocate(ffhh)
   deallocate(hice)
   deallocate(fice)
   deallocate(tisfc)
   deallocate(tprcp)
   deallocate(srflag)
   deallocate(snwdph)
   deallocate(slc)
   deallocate(shdmin)
   deallocate(shdmax)
   deallocate(slope)
   deallocate(snoalb)
   deallocate(oro,oro_uf)
   deallocate(sheleg)
   deallocate(sncovr)
   deallocate(hprime)
   deallocate(ozplin)
   deallocate(pl_lat,pl_pres,pl_time)
   deallocate( &
     SFALB   , &
     SFCDSW  , &
     COSZEN  , &
     PWAT    , &
     TMPMIN  , &
     TMPMAX  , &
     SPFHMIN , &
     SPFHMAX , &
     DUSFC   , &
     DVSFC   , &
     DTSFC   , &
     DQSFC   , &
     DLWSFC  , &
     ULWSFC  , &
     GFLUX   , &
     RUNOFF  , &
     EP      , &
     CLDWRK  , &
     DUGWD   , &
     DVGWD   , &
     PSMEAN  , &
     GESHEM  , &
     BENGSH  , &
     SFCNSW  , &
     SFCDLW  , &
     TSFLW   , &
     PSURF   , &
     U10M    , &
     V10M    , &
     HPBL    , &
     CHH     , &
     CMM     , &
     EPI     , &
     DLWSFCI , &
     ULWSFCI , &
     USWSFCI , &
     DSWSFCI , &
     DTSFCI  , &
     DQSFCI  , &
     GFLUXI  , &
     SRUNOFF , &
     T1      , &
     Q1      , &
     U1      , &
     V1      , &
     ZLVL    , &
     EVBSA   , &
     EVCWA   , &
     TRANSA  , &
     SBSNOA  , &
     SNOWCA  , &
     SOILM   , &
     SNOHFA  , &
     SMCWLT2 , &
     SMCREF2 , &
     suntim  , &                
     sfcemis , &                
     gsoil   , &
     gtmp2m  , &
     gustar  , &
     gpblh   , &
     gu10m   , &
     gv10m   , &
     gzorl   , &
     goro    )
  deallocate(cldcov,phy_f3d,phy_f2d)
  deallocate(ozjindx1,ozjindx2,ozddy)
  deallocate(hlw,swh,fluxr)
  deallocate(bfilt)
  is_allocated = .false.
 end subroutine destroy_phydata

 SUBROUTINE ORORD(LUGB,IORO,JORO,ORO,FNOROG)
 implicit none
 integer lugb, ioro, joro, kpdoro, ior, jor
 CHARACTER*(*) FNOROG
 real(r_single) oro4(ioro,joro), blnm, bltm
 real(r_kind) oro(ioro,joro)
 logical gausm
 kpdoro = 8
 IOR    = IORO
 JOR    = JORO
 ! fixrdg in sfcsub.f
 CALL FIXRDG(LUGB,IOR,JOR,FNOROG, &
             KPDORO,ORO4,GAUSM,BLNM,BLTM)
 if (ior .ne. ioro .or. jor .ne. joro) then
    print *,' orography file not o.k. run aborted'
    stop
 endif
 oro = oro4
 print *,'min/max oro:',minval(oro),maxval(oro)
 END SUBROUTINE ORORD

 SUBROUTINE FIXRDG(LUGB,IDIM,JDIM,FNGRIB, &
                  KPDS5,GDATA,GAUS,BLNO,BLTO)
 implicit none
 integer lgrib,n,lskip,jret,j,ndata,lugi,jdim,idim,lugb,&
         iret,kpds5,mdata,kdata
 CHARACTER*(*) FNGRIB
 PARAMETER(MDATA=2640*1320)
 REAL(r_single) GDATA(IDIM*JDIM)
 LOGICAL GAUS
 REAL(r_single) BLNO,BLTO
 real(r_single) data4(idim*jdim)
 LOGICAL*1 LBMS(MDATA)
 INTEGER KPDS(200),KGDS(200)
 INTEGER JPDS(200),JGDS(200), KPDS0(200)
 CLOSE(LUGB)
 call baopenr(lugb,fngrib,iret)
 IF (IRET .NE. 0) THEN
   WRITE(6,*) ' ERROR IN OPENING FILE ',trim(FNGRIB)
   PRINT *,'ERROR IN OPENING FILE ',trim(FNGRIB)
   CALL ABORT
 ENDIF
 WRITE(6,*) ' FILE ',trim(FNGRIB),&
 ' opened. Unit=',LUGB
 lugi    = 0
 lskip   = -1
 N       = 0
 JPDS    = -1
 JGDS    = -1
 JPDS(5) = KPDS5
 KPDS    = JPDS
 call getgbh(lugb,lugi,lskip,jpds,jgds,lgrib,ndata,&
             lskip,kpds,kgds,iret)
 WRITE(6,*) ' First grib record.'
 WRITE(6,*) ' KPDS( 1-10)=',(KPDS(J),J= 1,10)
 WRITE(6,*) ' KPDS(11-20)=',(KPDS(J),J=11,20)
 WRITE(6,*) ' KPDS(21-  )=',(KPDS(J),J=21,22)
 KPDS0=JPDS
 KPDS0(4)=-1
 KPDS0(18)=-1
 IF(IRET.NE.0) THEN
   WRITE(6,*) ' Error in GETGBH. IRET: ', iret
   IF (IRET == 99) WRITE(6,*) ' Field not found.'
   CALL ABORT
 ENDIF
 jpds = kpds0
 lskip = -1
 kdata=idim*jdim
 call getgb(lugb,lugi,kdata,lskip,jpds,jgds,ndata,lskip,&
                 kpds,kgds,lbms,data4,jret)
 if(jret.eq.0) then
   IF(NDATA.EQ.0) THEN
     WRITE(6,*) ' Error in getgb'
     WRITE(6,*) ' KPDS=',KPDS
     WRITE(6,*) ' KGDS=',KGDS
     STOP
   ENDIF
   IDIM=KGDS(2)
   JDIM=KGDS(3)
   gaus=kgds(1).eq.4
   blno=kgds(5)*1.e-3
   blto=kgds(4)*1.e-3
   !print *,'data min/max',minval(data4(1:idim*jdim)),maxval(data4(1:idim*jdim))
   gdata(1:idim*jdim)=data4(1:idim*jdim)
   WRITE(6,*) 'IDIM,JDIM=',IDIM,JDIM,&
   ' gaus=',gaus,' blno=',blno,' blto=',blto
 ELSE
   WRITE(6,*) ' Error in GETGB : JRET=',JRET
   WRITE(6,*) ' KPDS(13)=',KPDS(13),' KPDS(15)=',KPDS(15)
   STOP
 ENDIF
 END SUBROUTINE FIXRDG

 subroutine wrtout_sfc(fhour,filename)
   implicit none
   real(r_kind), intent(in) :: fhour
   character(len=120) filename
   integer,parameter :: version=200501
   type(sfcio_data) data
   integer iret,lu
   lu = 7

   ! sfchead is saved from input file as a private module variable.
   call sfcio_aldata(sfchead,data,iret)
   if(iret.ne.0) then
     write(6,*) 'error allocating surface data structure'
     stop
   endif
   sfchead%fhour   = fhour
   sfchead%idate   = idate_start
   data%tsea = tsea
   data%smc = smc
   data%sheleg = sheleg
   data%stc = stc
   data%tg3 = tg3
   data%zorl = zorl
   data%alvsf = alvsf
   data%alvwf = alvwf
   data%alnsf = alnsf
   data%alnwf = alnwf
   data%slmsk = slmsk
   data%vfrac = vfrac
   data%canopy = canopy
   data%f10m = f10m
   data%t2m = t2m
   data%q2m = q2m
   data%vtype = vtype
   data%stype = stype
   data%facsf = facsf
   data%facwf = facwf
   data%uustar = uustar
   data%ffmm = ffmm
   data%ffhh = ffhh
   data%hice = hice
   data%fice = fice
   data%tisfc = tisfc
!the addition of 8 Noah-related records starts here ...............
   data%tprcp = tprcp
   data%srflag = srflag
   data%snwdph = snwdph
   data%slc = slc
   data%shdmin = shdmin
   data%shdmax = shdmax
   data%slope = slope
   data%snoalb = snoalb
   data%orog = oro
! Not needed for version 200501
   !data%cv=cv
   !data%cvb=cvb
   !data%cvt=cvt
   call sfcio_swohdc(lu,filename,sfchead,data,iret)
   if(iret.ne.0) then
     write(6,*) 'error writing ',trim(filename)
     stop
   endif
   call sfcio_axdata(data,iret)

 end subroutine wrtout_sfc

 subroutine wrtout_flx(fhour,ta,filename)
   implicit none
   real(r_kind),intent(in) :: fhour
   real(r_double),intent(in) :: ta
   character(len=120),intent(in) :: filename
   integer, parameter :: noflx=7 
   real(r_kind) secswr,seclwr,zhour,fha,dtsw,dtlw
   integer,PARAMETER :: NFLD=25
   integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
   real(4) RTIMER(NFLD),rtime,rtimsw,rtimlw
   real(4) colat1
   integer,PARAMETER :: IPRS=1,ITEMP=11,IZNLW=33,IMERW=34,ISPHUM=51,IPWAT=54,&
            IPCPR=59,ISNOWD=65,ICLDF=71,ICCLDF=72,&
            ISLMSK=81,IZORL=83,IALBDO=84,ISOILM=144,ICEMSK=91,&
            ISIK=92,                               &
            ILHFLX=121,ISHFLX=122,IZWS=124,IMWS=125,IGHFLX=155,&
            IUSWFC=160,IDSWFC=161,IULWFC=162,IDLWFC=163,&
            INSWFC=164,INLWFC=165,&
            IDSWVB=166,IDSWVD=167,IDSWNB=168,IDSWND=169,&
            ITMX=15,ITMN=16,IRNOF=90,IEP=145,&
            IQMX=204,IQMN=205,&
            ICLDWK=146,IZGW=147,IMGW=148,IHPBL=221,&
            IDSWF=204,IDLWF=205,IUSWF=211,IULWF=212,ICPCPR=214,&
            IUVBF=200,IUVBFC=201,ISUNTM=191,&
            ICSUSW=160,ICSDSW=161,ICSULW=162,ICSDLW=163,&
            ISFC=1,ITOA=8,IELEV=105,&
            ISGLEV=109,IDBLS=111,I2DBLS=112,ICOLMN=200,&
            IBLBL=209,IBLTL=210,IBLLYR=211,&
            ILCBL=212,ILCTL=213,ILCLYR=214,&
            IMCBL=222,IMCTL=223,IMCLYR=224,&
            IHCBL=232,IHCTL=233,IHCLYR=234,&
            ICVBL=242,ICVTL=243,ICVLYR=244,&
            ISLC=160,ISNOD=66,&
            ISLO=222,ISBS=198,ISNC=238,ICMM=179,&
            ISNOHF=229,ISMCWLT=219,ISMCREF=220,&
            ICNP=223,&
            IVEG=87,IVTP=225,ISTP=224,IUST=253,IHGT=7,&
            IRST=140,ICHH=208,ISRF=235,IEVBS=199,&
            IEVCW=200,ITRAN=210,ISTC=86,&
            INST=10,IWIN=2,IAVG=3,IACC=4,&
            ICEN=7,IFHOUR=1,IFDAY=2
   LOGICAL(1) LBM(nlons*nlats)
   CHARACTER G(200+nlons*nlats*(16+1)/8)
   INTEGER   IPUR(NFLD),ITLR(NFLD),ICEN2,IGEN
   DATA      IPUR/IULWF , IUSWF , IUSWF , IDSWF ,  ICLDF,   IPRS,&
                  IPRS, ITEMP ,  ICLDF,   IPRS,   IPRS, ITEMP ,&
                  ICLDF,   IPRS,   IPRS, ITEMP ,  IUVBF, IUVBFC ,&
                  IDSWF, ICSULW, ICSUSW, ICSDLW, ICSUSW, ICSDSW,&
                  ICSULW /
   DATA      ITLR/ITOA  , ITOA  , ISFC  , ISFC  , IHCLYR, IHCTL ,&
                  IHCBL , IHCTL , IMCLYR, IMCTL , IMCBL , IMCTL ,&
                  ILCLYR, ILCTL , ILCBL , ILCTL , ISFC  , ISFC ,&
                  ITOA  ,  ITOA ,  ITOA ,  ISFC , ISFC  , ISFC,&
                  ISFC /
   INTEGER   IENS(5),IDS(255),i,j,n,k,k4,l
   real(4) wrkga(nlons*nlats)
   integer slmski(nlons*nlats)

   call BAOPENWT(NOFLX,filename,ierr)
   if (ierr .ne. 0) print *,' iostat after baclose of flxf ile ',ierr

   fha = ta/3600.
   zhour = fhour - fha ! last time accum arrays zeroed
   !print *,'fha,fhour,zhour=',fha,fhour,zhour
   if (fhswr .lt. dt) then
      dtsw = dt
   else
      dtsw  = 3600.0 * fhswr
   endif
   if (fhlwr .lt. dt) then
      dtlw = dt
   else
      dtlw  = 3600.0 * fhlwr
   endif
   SECSWR=MAX(ta,DTSW)
   SECLWR=MAX(ta,DTLW)
   ICEN2 = sighead%icen2; IGEN = sighead%igen

   CALL IDSDEF(1,IDS)
   ids(IQMX)   = 5
   ids(IQMN)   = 5
   ids(IUVBF)  = 2
   ids(IUVBFC) = 2
   ids(icemsk) = 3      ! ICE CONCENTRATION ()
   ids(isik)   = 2      ! ICE THICKNESS (M)
   ids(IZORL)  = 4
   ids(IHGT)   = 3
   ids(IVEG)   = 2
   ids(IUST)   = 3
   ids(ICHH)   = 4
   ids(ICMM)   = 4
   ids(ISRF)   = 5
   ids(ITEMP)  = 3
   ids(ISPHUM) = 6
   ids(IZNLW)  = 2
   ids(IMERW)  = 2
   ids(ISNC)   = 3
   ids(ISTC)   = 4
   ids(ISOILM) = 4
   ids(ISNOD)  = 6
   ids(ISNOWD) = 5
   ids(ICNP)   = 5
   ids(IPCPR)  = 6
   ids(ICPCPR) = 6
   ids(IRNOF)  = 5                                            
   ids(ISMCWLT)  = 4
   ids(ISMCREF)  = 4
   ILPDS   = 28
   IF(ICEN2.EQ.2) ILPDS=45
   IENS(1) = 1
   IENS(2) = sighead%iens(1)
   IENS(3) = sighead%iens(2)
   IENS(4) = 1
   IENS(5) = 255
   IYR     = IDATE_START(4)
   IMO     = IDATE_START(2)
   IDA     = IDATE_START(3)
   IHR     = IDATE_START(1)
   IFHR    = NINT(ZHOUR)
   ITHR     = NINT(FHOUR)
   IF(FHOUR.GT.ZHOUR) THEN
     RTIME = 1./(3600.*(FHOUR-ZHOUR))
   ELSE
     RTIME = 0.
   ENDIF
   IF(SECSWR.GT.0.) THEN
     RTIMSW = 1./SECSWR
   ELSE
     RTIMSW = 1.
   ENDIF
   IF(SECLWR.GT.0.) THEN
     RTIMLW = 1./SECLWR
   ELSE
     RTIMLW = 1.
   ENDIF
   RTIMER    = RTIMSW
   RTIMER(1) = RTIMLW
   RTIMER(20)=RTIMLW       ! CSULF_TOA
   RTIMER(22)=RTIMLW       ! CSDLF_SFC
   RTIMER(25)=RTIMLW       ! CSULF_SFC

   ! FIRST COLATITUDE OF GRID IF IDRT=4 (RADIANS)
   colat1 = 0.5*con_pi-lats(1,1)
    
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      slmski(n) = nint(slmsk(i,j))
   enddo
   LBM = slmski.EQ.1

   call twodtooned(dusfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IZWS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IZWS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '01)Zonal compt of momentum flux (N/m**2) land and sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dvsfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IMWS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IMWS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '02)Merid compt of momentum flux (N/m**2) land and sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dtsfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '03)Sensible heat flux (W/m**2) land and sea surface           '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dqsfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ILHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ILHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
   '04)Latent heat flux (W/m**2) land and sea surface             '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(tsea,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IMWS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '05)Temperature (K) land and sea surface                       '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(smc(:,:,1),wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISOILM,I2DBLS,0,10,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '06)Volumetric soil moist content (frac) layer 10cm and 0cm    '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(smc(:,:,2),wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISOILM,I2DBLS,10,40,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
   '07)Volumetric soil moist content (frac) layer 40cm and 10cm  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(stc(:,:,1),wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ITEMP,I2DBLS,0,10,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
   '08)Temp (K) layer betw two depth below land sfc 10cm and 0cm  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(stc(:,:,2),wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ITEMP,I2DBLS,10,40,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
   '09)Temp (K) layer betw two depth below land sfc 40cm and 10cm'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(sheleg,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISNOWD,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISNOWD),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '10)Water equiv of accum snow depth (kg/m**2) land sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dlwsfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IDLWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IDLWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '11)Downward long wave radiation flux (W/m**2) land sea surface'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(ulwsfc*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IULWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IULWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '12)Upward long wave radiation flux (W/m**2) land sea surface  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
!..........................................................
!.......  FIX FLUXES FOR APPROX DIURNAL CYCLE
   do k=1,4
      call twodtooned(fluxr(:,:,k)*rtimer(k),wrkga)
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                  0,IPUR(K),ITLR(K),0,0,IYR,IMO,IDA,IHR,&
                  IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(K)),IENS,&
                  0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.1)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '13)Upward long wave radiation flux (W/m**2) top of atmosphere '
      if(ierr.ne.0.and.k.eq.2)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '14)Upward solar radiation flux (W/m**2) top of atmosphere     '
      if(ierr.ne.0.and.k.eq.3)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '15)Upward solar radiation flux (W/m**2) land and sea surface  '
      if(ierr.ne.0.and.k.eq.4)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '16)Downward solar radiation flux (W/m**2) land and sea surface'
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
   enddo       
!..........................................................
!
!     For UV-B fluxes
!
   call twodtooned(fluxr(:,:,21)*rtimsw,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,129,ICEN,IGEN,&
               0,IUVBF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IUVBF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '17)UV-B Downward solar flux (W/m**2) land sea surface'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(fluxr(:,:,22)*rtimsw,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,129,ICEN,IGEN,&
               0,IUVBFC,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IUVBFC),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '18)clear sky UV-B Downward solar flux (W/m**2) land sea surface'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   do k=5,7
      call twodtooned(fluxr(:,:,k)*100.*rtimsw,wrkga)
      K4=4+(K-5)*4
      L=K4+1
      LBM=wrkga.ge.0.5
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                   0,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,&
                   IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,&
                   0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '19)Total cloud cover (percent) high cloud layer               '
      if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '23)Total cloud cover (percent) middle cloud layer             '
      if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '27)Total cloud cover (percent) low cloud layer                '
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
 
      wrkga = 0.
      do n=1,nlons*nlats
         j = 1+(n-1)/nlons
         i = n-(j-1)*nlons
         if (fluxr(i,j,k) .gt. 0) then
            wrkga(n) = fluxr(i,j,k+3)*1000./fluxr(i,j,k)
         endif
      enddo
      L=K4+2
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                   1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,&
                   IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,&
                   0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '20)Pressure (Pa) high cloud top level                         '
      if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '24)Pressure (Pa) middle cloud top level                       '
      if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '28)Pressure (Pa) low cloud top level                          '
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
 
      wrkga = 0.
      do n=1,nlons*nlats
         j = 1+(n-1)/nlons
         i = n-(j-1)*nlons
         if (fluxr(i,j,k) .gt. 0) then
            wrkga(n) = fluxr(i,j,k+6)*1000./fluxr(i,j,k)
         endif
      enddo
      L=K4+3
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                   1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,&
                   IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,&
                   0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '21)Pressure (Pa) high cloud bottom level                      '
      if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '25)Pressure (Pa) middle cloud bottom level                    '
      if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '29)Pressure (Pa) low cloud bottom level                       '
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
 
      wrkga = 0.
      do n=1,nlons*nlats
         j = 1+(n-1)/nlons
         i = n-(j-1)*nlons
         if (fluxr(i,j,k) .gt. 0) then
            wrkga(n) = fluxr(i,j,k+9)*1000./fluxr(i,j,k)
         endif
      enddo
      L=K4+4
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                   1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,&
                   IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,&
                   0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '22)Temperature (K) high cloud top level                       '
      if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '26)Temperature (K) middle cloud top level                     '
      if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',&
       '30)Temperature (K) low cloud top level                        '
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
   enddo

   call twodtooned(1000.*geshem*rtime,wrkga)
   !print *,'min/max prate',minval(wrkga),maxval(wrkga),rtime
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IPCPR,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPCPR),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '31)Precipitation rate (kg/m**2/s) land and sea surface        '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(1000.*bengsh*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICPCPR,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICPCPR),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '32)Convective precipitation rate (kg/m**2/s) land sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   LBM = slmski.EQ.1
   call twodtooned(gflux*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IGHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IGHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '33)Ground heat flux (W/m**2) land and sea surface             '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(MOD(slmsk,2.),wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISLMSK,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISLMSK),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '34)Land-sea mask (1=land; 0=sea) (integer) land sea surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(fice,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICEMSK,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICEMSK),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '35)Ice concentration (ice>0; no ice=0) (1/0) land sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(u10m,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
              0,IZNLW,IELEV,0,10,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZNLW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '36)u wind (m/s) height above ground                           '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(v10m,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IMERW,IELEV,0,10,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IMERW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '37)v wind (m/s) height above ground                           '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(t2m,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ITEMP,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '38)Temperature (K) height above ground                        '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(q2m,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISPHUM,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISPHUM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '39)Specific humidity (kg/kg) height above ground              '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(psurf,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IPRS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '40)Pressure (Pa) land and sea surface                         '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(tmpmax,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ITMX,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(ITMX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '41)Maximum temperature (K) height above ground                '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(tmpmin,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ITMN,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(ITMN),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '42)Minimum temperature (K) height above ground                '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(spfhmax,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,133,ICEN,IGEN,&
               0,IQMX,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(IQMX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '41a)Maximum specific humidity (kg/kg) height above ground      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(spfhmin,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,133,ICEN,IGEN,&
               0,IQMN,IELEV,0,2,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(IQMN),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '42a)Minimum specific humidity (kg/kg) height above ground      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(1000.*runoff,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IRNOF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IACC,0,0,ICEN2,IDS(IRNOF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '43)Runoff (kg/m**2) land and sea surface                      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(ep*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IEP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '44)Potential evaporation rate (w/m**/) land and sea surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(cldwrk*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICLDWK,ICOLMN,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDWK),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '45)Cloud work function (J/Kg) total atmospheric column        '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dugwd*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IZGW,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IZGW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '46)Zonal gravity wave stress (N/m**2) land and sea surface    '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dvgwd*rtime,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IMGW,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IMGW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '47)Meridional gravity wave stress (N/m**2) land sea surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(hpbl,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IHPBL,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IHPBL),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '48)Boundary layer height '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(pwat,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IPWAT,ICOLMN,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPWAT),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '49)Precipitable water (kg/m**2) total atmospheric column      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   wrkga = 0.
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      if (fluxr(i,j,4) .gt. 0) then
         wrkga(n) = fluxr(i,j,3)/fluxr(i,j,4) * 100.
         if (wrkga(n) .gt. 100.) wrkga(n) = 100.
      endif
   enddo
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IALBDO,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IALBDO),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '50)Albedo (percent) land and sea surface                      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(fluxr(:,:,26)*100.*rtimsw,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICLDF,ICOLMN,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '51)Total cloud cover (percent) total atmospheric column       '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
!
! CONVECTIVE CLOUDS
! LABELED INSTANTANEOUS BUT ACTUALLY AVERAGED OVER FHSWR HOURS
!
   call twodtooned(cv*1.e12,wrkga)
   LBM=wrkga.Ge.0.5
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICLDF,ICVLYR,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICLDF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '52)Total cloud cover (percent) convective cloud layer         '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

!  wrkga = 0.
!  do n=1,nlons*nlats
!     j = 1+(n-1)/nlons
!     i = n-(j-1)*nlons
!     if (cv(i,j) .gt. 0) then
!        wrkga(n) = cvt(i,j)*1.e3
!     endif
!  enddo
!  call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
!              1,IPRS,ICVTL,0,0,IYR,IMO,IDA,IHR,&
!              IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,&
!              0.,0.,0.,0.,0.,0.,G,LG,IERR)
!  if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
!   '53)Pressure (Pa) convective cloud top level                   '
!  IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

!  wrkga = 0.
!  do n=1,nlons*nlats
!     j = 1+(n-1)/nlons
!     i = n-(j-1)*nlons
!     if (cv(i,j) .gt. 0) then
!        wrkga(n) = cvb(i,j)*1.e3
!     endif
!  enddo
!  call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
!              1,IPRS,ICVBL,0,0,IYR,IMO,IDA,IHR,&
!              IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,&
!              0.,0.,0.,0.,0.,0.,G,LG,IERR)
!  if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
!   '54)Pressure (Pa) convective cloud bottom level                '
!  IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

!.................................................
!...   SAVE B.L. CLOUD AMOUNT
!
   call twodtooned(fluxr(:,:,27)*100.*rtimsw,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                 0,ICLDF,IBLLYR,0,0,IYR,IMO,IDA,IHR,&
                 IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDF),IENS,&
                 0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '55)Total cloud cover (percent) boundary layer cloud layer     '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(hice,wrkga)
   LBM=slmski.NE.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISIK,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISIK),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '56)Sea ice thickness (m) category 1'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(smc(:,:,3),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                1,ISOILM,I2DBLS,40,100,IYR,IMO,IDA,IHR,&
                IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
                0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '57)Volumetric soil moist content (frac) layer 100cm and 40cm '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(smc(:,:,4),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISOILM,I2DBLS,100,200,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '58)Volumetric soil moist content (frac) layer 200cm and 100cm '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(stc(:,:,3),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ITEMP,I2DBLS,40,100,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '59)Temp (K) layer betw two depth below land sfc 100cm and 40cm'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(stc(:,:,4),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ITEMP,I2DBLS,100,200,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '60)Temp (K) layer betw two depth below land sfc 200cm and 100cm'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(slc(:,:,1),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISLC,I2DBLS,0,10,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '61)Liquid soil moist content (frac) layer 10cm and 0cm  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(slc(:,:,2),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISLC,I2DBLS,10,40,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '62)Liquid soil moist content (frac) layer 40cm and 10cm '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(slc(:,:,3),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISLC,I2DBLS,40,100,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '63)Liquid soil moist content (frac) layer 100cm and 40cm'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(slc(:,:,4),wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISLC,I2DBLS,100,200,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
   '64)Liquid soil moist content (frac) layer 200cm and 100cm'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(snwdph/1.e3,wrkga)
   LBM=slmski.EQ.1 .or. slmski.EQ.2
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISNOD,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISNOD),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '65)Snow depth (m) land surface                  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(canopy,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ICNP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICNP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '66)Canopy water content (kg/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(zorl/100.,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
              0,IZORL,ISFC,0,0,IYR,IMO,IDA,IHR,&
              IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZORL),IENS,&
              0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '67)Surface roughness (m)    '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(vfrac*100.,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
              1,IVEG,ISFC,0,0,IYR,IMO,IDA,IHR,&
              IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IVEG),IENS,&
              0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '68)Vegetation fraction (fractional) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(vtype,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IVTP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IVTP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '69)Vegetation type land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(stype,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISTP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISTP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '70)Soil type land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(slope,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISLO,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISLO),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '71)Slope type land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(uustar,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IUST,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IUST),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '72)Frictional velocity (m/s)'
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(oro,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IZWS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IZWS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '73)Surface height (m)       '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(srflag,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
              1,IRST,ISFC,0,0,IYR,IMO,IDA,IHR,&
              IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IRST),IENS,&
              0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '74)Freezing precip flag land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(chh,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ICHH,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICHH),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ', &
    '75)Exchange coefficient CH(m/s)       '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(cmm,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               0,ICMM,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICMM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '76)Exchange coefficient CM(m/s)      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(epi,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IEP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IEP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '77)Potential evaporation rate (w/m**2) land and sea surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dlwsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IDLWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IDLWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '78)Downward long wave radiation flux (W/m**2) '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(ulwsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IULWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IULWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '79)Upward long wave radiation flux (W/m**2)  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(uswsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IUSWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IUSWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '80)Upward short wave radiation flux (W/m**2)  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dswsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IDSWF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IDSWF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '81)Downward short wave radiation flux (W/m**2)  '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(dtsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   IF(IERR.EQ.0) then
      CALL WRYTE(noflx,LG,G)
   else
       print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '82)Sensible heat flux (W/m**2) land and sea surface          '
   endif

   call twodtooned(dqsfci,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ILHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ILHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   IF(IERR.EQ.0) then
     CALL WRYTE(noflx,LG,G)
   else
     print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '83)Latent heat flux (W/m**2) land and sea surface            '
   endif

   call twodtooned(gfluxi,wrkga)
   LBM=slmski.NE.0
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IGHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IGHFLX),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   IF(IERR.EQ.0) then
     CALL WRYTE(noflx,LG,G)
   else
     print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '84)Ground heat flux (W/m**2) land and sea surface            '
   endif

   call twodtooned(srunoff*1000.,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISRF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISRF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '85)Surface runoff (kg/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(t1,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ITEMP,isglev,1,1,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '86)Lowest model level Temp (K)       '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(q1,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,ISPHUM,isglev,1,1,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISPHUM),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '87)Lowest model specific humidity (kg/kg)      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(u1,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IZNLW,isglev,1,1,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZNLW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '88)Lowest model u wind (m/s)      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(v1,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IMERW,isglev,1,1,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IMERW),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '89)Lowest model v wind (m/s)      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(zlvl,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IHGT,isglev,1,1,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IHGT),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '90)Lowest model level height (m) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(evbsa*rtime,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IEVBS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEVBS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '91)Direct evaporation from bare soil(W/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(evcwa*rtime,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,IEVCW,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEVBS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '92)Canopy water evaporation(W/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(transa*rtime,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ITRAN,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ITRAN),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '93)Transpiration (W/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(sbsnoa*rtime,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISBS,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISBS),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '94)Snow Sublimation (W/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(snowca*rtime*100.,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISNC,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISNC),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '95)Snow Cover (fraction) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(soilm*1000.,wrkga) !! convert from m to (mm)kg/m^2
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               1,ISTC,I2DBLS,0,200,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISTC),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '96)Total column soil moisture (Kg/m^2) land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   do k=19,25
      if (k .eq. 19) then
        L = 18
      else
        L = k + 8
      endif
      call twodtooned(fluxr(:,:,L)*rtimer(k),wrkga)
      call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
                  0,IPUR(K),ITLR(K),0,0,IYR,IMO,IDA,IHR,&
                  IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(K)),IENS,&
                  0.,0.,0.,0.,0.,0.,G,LG,IERR)
      if(ierr.ne.0.and.k.eq.19)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '97)Downward solar radiation flux (W/m**2) TOA '
      if(ierr.ne.0.and.k.eq.20)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '98)CS upward long wave radiation flux (W/m**2) TOA '
      if(ierr.ne.0.and.k.eq.21)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '99)CS upward solar radiation flux (W/m**2) TOA     '
      if(ierr.ne.0.and.k.eq.22)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '100)CS downward long radiation flux (W/m**2) SFC  '
      if(ierr.ne.0.and.k.eq.23)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '101)CS upward solar radiation flux (W/m**2)  SFC '
      if(ierr.ne.0.and.k.eq.24)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '102)CS downward solar radiation flux (W/m**2) SFC'
      if(ierr.ne.0.and.k.eq.25)print*,'wrtsfc gribit ierr=',ierr,'  ',&
      '103)CS upward long wave radiation flux (W/m**2) SFC '
      IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
   enddo

   call twodtooned(snohfa*rtime,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISNOHF,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISNOHF),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '104)Snow phase-change heat flux [W/m^2] land surface      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(smcwlt2,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
               1,ISMCWLT,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISMCWLT),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   IF(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '105)Wilting point [fraction] land surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
!..........................................................
   call twodtooned(smcref2,wrkga)
   LBM=slmski.EQ.1
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,130,ICEN,IGEN,&
             1,ISMCREF,ISFC,0,0,IYR,IMO,IDA,IHR,&
             IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISMCREF),IENS,&
             0.,0.,0.,0.,0.,0.,G,LG,IERR)
   IF(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '106)Field capacity [fraction] land surface   '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)
 
!  call twodtooned(suntim,wrkga)
!  call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,133,ICEN,IGEN,&
!              0,ISUNTM,ISFC,0,0,IYR,IMO,IDA,IHR,&
!              IFHOUR,IFHR,ITHR,IACC,0,0,ICEN2,IDS(ISUNTM),IENS,&
!              0.,0.,0.,0.,0.,0.,G,LG,IERR)
!  if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
!   '107)Accumulated sunshine duration time (sec) '
!  IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   PRINT *,'GRIB FLUX FILE WRITTEN ',fhour,idate_start
   call baclose(noflx,ierr)
   print *,' iostat after baclose of flxf file ',ierr

 end subroutine wrtout_flx

 subroutine twodtooned(data2,data1)
   real(r_kind), intent(in) :: data2(nlons,nlats)
   real(4), intent(out) :: data1(nlons*nlats)
   integer i,j,n
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      data1(n) = data2(i,j)
   enddo
 end subroutine twodtooned

end module phy_data
