 module phy_run
! compute physics tendencies for gfs physics
! getphytend: compute total physics tendencies.
! Public subroutines:
! getphytend: compute tendencies in spectral space.

 use params, only: nlevs,nlons,nlats,ntrunc,ndimspec,ntrac,nmtvr,idate_start,dt,&
 fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,&
 ncw,iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,&
 psautco,mstrat,pre_rad,bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,ntoz,ntclw
 use kinds, only: r_kind
 use shtns, only: grdtospec, getvrtdivspec, lons, lats
 use grid_data, only: virtempg,dlnpdtg,tracerg,ug,vg
 use pressure_data, only:  prs,psg,pk,ak,bk
 use phy_data, only: solcon,slag,sdec,cdec,nfxr,ncld,&
    lsoil,timeoz,latsozp,levozp,pl_coeff,ozplin,pl_pres,pl_time,&
    slmsk,sheleg,sncovr,snoalb,zorl,hprime,alvsf,ozjindx1,ozjindx2,ozddy,&
    alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,coszen,cv,cvt,cvb,sfcemis,&
    sfcnsw,sfcdsw,sfalb,sfcdlw,tsflw,slc,smc,stc,slope,shdmin,shdmax, &
    vfrac,tg3,stype,vtype,oro,oro_uf,uustar,phy_f3d,phy_f2d,&
    fluxr,hlw,swh, &
    hice  ,    fice  ,        &
    tisfc ,    tsea  ,        &
    tprcp ,    cv    ,        &
    cvb   ,    cvt   ,        &
    srflag,    snwdph,        &
    sheleg,    sncovr,        &
    zorl  ,    canopy,        &
    ffmm  ,    ffhh  ,        &
    f10m  ,    srunoff,       &
    evbsa ,    evcwa ,        &
    snohfa,    transa,        &
    sbsnoa,    snowca,        &
    soilm ,    tmpmin,        &
    tmpmax,    dusfc ,        &
    dvsfc ,    dtsfc ,        &
    dqsfc ,    geshem,        &
    gflux ,    dlwsfc,        & 
    ulwsfc,    suntim,        &
    runoff,    ep    ,        &
    cldwrk,    dugwd ,        &
    dvgwd ,    psmean,        &
    bengsh,    spfhmin,  spfhmax,      &
    t2m   ,    q2m   ,        &
    u10m  ,    v10m  ,        &
    zlvl  ,    psurf ,        &
    hpbl  ,    pwat  ,        &
    t1    ,    q1    ,        &
    u1    ,    v1    ,        &
    chh   ,    cmm   ,        &
    dlwsfci,   ulwsfci,       &
    dswsfci,   uswsfci,       &
    dtsfci,    dqsfci,        &
    gfluxi,    epi   ,        &
    smcwlt2,   smcref2,       &
    gsoil,     gtmp2m,        &
    gustar,    gpblh,         &
    gu10m,     gv10m,         &
    gzorl,     goro          
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
               eps => con_eps, omega => con_omega, cvap => con_cvap, &
               grav => con_g

 implicit none
 private

 public :: getphytend

 contains

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,t,dtx)
   use omp_lib, only: omp_get_num_threads, omp_get_thread_num
   use mersenne_twister, only : random_setseed, random_index, random_stat
   use module_radiation_astronomy, only: astronomy
   use module_radiation_driver, only: radinit, grrad
   ! compute physics tendencies for gfs physics
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
   dtracerspecdt
   real(r_kind), intent(in) :: t,dtx
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   real(r_kind), parameter :: qmin = 1.e-10 ! min value for clipping tracers
   real(r_kind), parameter :: typical_pgr = 95000.0
   real(r_kind)  :: fscav(ntrac-ncld-1),ozp(levozp,pl_coeff),&
   fhour,dtsw,dtlw,facoz,clstp,solhr,dphi,dpshc(1),delta,rk,&
   gt(nlevs),prsl(nlevs),prsi(nlevs+1),vvel(nlevs),&
   f_ice(nlevs),f_rain(nlevs),r_rime(nlevs),&
   prslk(nlevs),gq(nlevs,ntrac),prsik(nlevs+1),si_loc(nlevs+1),phii(nlevs+1),&
   phil(nlevs),rann(1),acv(1),acvb(1),acvt(1),adt(nlevs),adu(nlevs),adv(nlevs),&
   dum1(1),rqtk(1),upd_mf(nlevs),dwn_mf(nlevs),det_mf(nlevs),dkh(nlevs),rnp(nlevs),&
   adq(nlevs,ntrac),dt3dt(nlevs,6),du3dt(nlevs,4),dv3dt(nlevs,4),dq3dt(nlevs,5+pl_coeff),&
   phy3d(nlevs,num_p3d),phy2d(num_p2d),hlw_tmp(nlevs),swh_tmp(nlevs),fluxr_tmp(nfxr),&
   cldcov_tmp(nlevs),hprime_tmp(nmtvr),slc_tmp(lsoil),smc_tmp(lsoil),stc_tmp(lsoil),&
   gu(nlevs),gv(nlevs)
   integer :: icsdsw(1),icsdlw(1),ilons(1)
   real(r_kind), allocatable, dimension(:,:) :: coszdg,dpsdt
   real(r_kind), allocatable, dimension(:,:,:) :: ozplout,dtdt,dudt,dvdt
   real(r_kind), allocatable, dimension(:,:,:) :: cldcov
   real(r_kind), allocatable, dimension(:,:,:,:) :: dtracersdt
   logical, save :: sas_shal
   logical, save :: ipsd0
   logical, save :: first=.true.
   logical :: lsswr=.false. ! sw rad call
   logical :: lslwr=.false. ! sw rad call
   logical :: lssav=.true. ! store 3d cloud field?
   logical :: lssav_cc=.false. ! flag for coupling
   logical :: flipv=.true. ! vert dir flip for RAS
   logical :: ldiag3d = .false.
   logical :: lggfs3d = .false.
   integer :: idat(8), jdat(8)
   type (random_stat) :: stat
   integer :: numrdm(nlons*nlats*2), ixseed(nlons,nlats,2)
   integer :: i,j,k,n,nt,ipseed,nstep,nswr,nlwr,nnrcm
   integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds
   real(8) tstart,tend
   integer(8) count, count_rate, count_max

   rk = rd/cp
   delta = cvap/cp-1. ! used for virt temp to sensible temp computation

   allocate(ozplout(levozp,nlats,pl_coeff))
   allocate(coszdg(nlons,nlats))
   allocate(cldcov(nlons,nlats,nlevs)) ! clouds diagnosed by grrad (diagnostic)
   allocate(dtdt(nlons,nlats,nlevs),dudt(nlons,nlats,nlevs),dvdt(nlons,nlats,nlevs))
   allocate(dtracersdt(nlons,nlats,nlevs,ntrac))
   allocate(dpsdt(nlons,nlats))

! is it a radiation time step (long or short wave)?
   fhour = t/3600.
   nswr = int(fhswr*3600./dt)
   nlwr = int(fhlwr*3600./dt)
   if (t > 0) then
      nstep = int(t/dt)
      lsswr = mod(nswr,nstep) .eq. 0
      lslwr = mod(nlwr,nstep) .eq. 0
   else
      nstep=0
      lsswr=.true.
      lslwr=.true.
   endif
! compute forecast valid time (jdat).
   call getvaliddate(fhour,idate_start,idat,jdat)

! initialize on first call
   if (first) then
     sas_shal = sashal .and. (.not. ras)
! get a sigma distribution for radiation-cloud initialization
     si_loc(nlevs+1)= ak(1)/typical_pgr+bk(1)
     do k=1,nlevs
        si_loc(nlevs+1-k)= ak(k+1)/typical_pgr + bk(k+1)
        !print *,k,si_loc(nlevs+1-k)
     enddo
! generate initial permutation seed for random number generator
     if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
        ipsd0 = 17*idate_start(1) + 43*idate_start(2) + &
                37*idate_start(3) + 23*idate_start(4)
        print *,'  Radiation sub-cloud initial seed =',ipsd0,&
                ' idate =',idate_start
     endif
     first = .false.
   endif         ! end_if_first
!
!===> *** ...  radiation initialization
!
   if (lsswr .or. lslwr) then ! call radiation
! the following block of code is adapted from gloopr.f
   dtsw  = 3600.0 * fhswr
   dtlw  = 3600.0 * fhlwr
   call radinit                                                      &
!  ---  input:
     &     ( si_loc, NLEVS, 1, idat, jdat, ICTM, ISOL, ICO2,         &
     &       IAER, IALB, IEMS, 1, NUM_P3D, ISUBC_SW, ISUBC_LW,       &
     &       IOVR_SW, IOVR_LW, 0 )
!  ---  output: ( none )
!===> *** ...  astronomy for sw radiation calculation.
!
  call astronomy                                                    &
!  ---  inputs:
       ( nlons, nlats, lons, lats,                                  &
         fhswr, jdat, lsswr,                                        &
!  ---  outputs:
         solcon, slag, sdec, cdec, coszen, coszdg                   &
        )
!
!===> *** ...  generate 2-d random seeds array for sub-grid cloud-radiation
!
      if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
        ipseed = mod(nint(100.0*sqrt(fhour*3600)), ipsdlim) + 1 + ipsd0
        call random_setseed                                             &
!  ---  inputs:
           ( ipseed,                                                    &
!  ---  outputs:
             stat                                                       &
            )
        call random_index                                               &
!  ---  inputs:
           ( ipsdlim,                                                   &
!  ---  outputs:
             numrdm, stat                                               &
           )
        do k = 1, 2
          do j = 1, nlats
            do i = 1, nlons
              ixseed(i,j,k) = numrdm(i+(j-1)*nlons+(k-1)*nlons*nlats)
            enddo
          enddo
        enddo
      endif

! compute radiation tendencies.
! control parameters (all set in namelist): 

!  ---  ICTM=yyyy#, controls time sensitive external data (e.g. CO2, solcon, aerosols, etc)
!     integer, parameter :: ICTM =   -2 ! same as ICTM=0, but add seasonal cycle from
!                                       ! climatology. no extrapolation.
!     integer, parameter :: ICTM =   -1 ! use user provided external data set for the
!                                       ! forecast time. no extrapolation.
!     integer, parameter :: ICTM =    0 ! use data at initial cond time, if not
!                                       ! available, use latest, no extrapolation.
!!    integer, parameter :: ICTM =    1 ! use data at the forecast time, if not
!                                       ! available, use latest and extrapolation.
!     integer, parameter :: ICTM =yyyy0 ! use yyyy data for the forecast time,
!                                       ! no further data extrapolation.
!     integer, parameter :: ICTM =yyyy1 ! use yyyy data for the fcst. if needed, do
!                                       ! extrapolation to match the fcst time.

!  ---  ISOL controls solar constant data source
!!    integer, parameter :: ISOL = 0   ! use prescribed solar constant
!     integer, parameter :: ISOL = 1   ! use varying solar const with 11-yr cycle

!  ---  ICO2 controls co2 data source for radiation
!     integer, parameter :: ICO2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ICO2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ICO2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  IALB controls surface albedo for sw radiation
!!    integer, parameter :: IALB = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: IALB = 1   ! use modis derived alb (to be developed)

!  ---  IEMS controls surface emissivity and sfc air/ground temp for lw radiation
!        ab: 2-digit control flags. a-for sfc temperature;  b-for emissivity
!!    integer, parameter :: IEMS = 00  ! same air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 01  ! same air/ground temp; varying veg typ based emis
!!    integer, parameter :: IEMS = 10  ! diff air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 11  ! diff air/ground temp; varying veg typ based emis

!  ---  IAER  controls aerosols scheme selections
!     Old definition
!     integer, parameter :: IAER  = 1  ! opac climatology, without volc forcing
!     integer, parameter :: IAER  =11  ! opac climatology, with volcanic forcing
!     integer, parameter :: IAER  = 2  ! gocart prognostic, without volc forcing
!     integer, parameter :: IAER  =12  ! gocart prognostic, with volcanic forcing
!     New definition in this code IAER = abc (a:volcanic; b:lw; c:sw)
!                             b, c values: (0:none; 1:opac; 2:gocart)
!  IAER =   0 --> no aerosol effect at all (volc, sw, lw)
!       =   1 --> only tropospheric sw aerosols, no trop-lw and volc
!       =  10 --> only tropospheric lw aerosols, no trop-sw and volc
!       =  11 --> both trop-sw and trop-lw aerosols, no volc
!       = 100 --> only strato-volc aeros, no trop-sw and trop-lw
!       = 101 --> only sw aeros (trop + volc), no lw aeros
!       = 110 --> only lw aeros (trop + volc), no sw aeros
!       = 111 --> both sw and lw aeros (trop + volc) 

!  ---  IOVR controls cloud overlapping method in radiation:
!     integer, parameter :: IOVR_SW = 0  ! sw: random overlap clouds
!!    integer, parameter :: IOVR_SW = 1  ! sw: max-random overlap clouds

!     integer, parameter :: IOVR_LW = 0  ! lw: random overlap clouds
!!    integer, parameter :: IOVR_LW = 1  ! lw: max-random overlap clouds

!  ---  ISUBC controls sub-column cloud approximation in radiation:
!     integer, parameter :: ISUBC_SW = 0 ! sw: without sub-col clds approx
!     integer, parameter :: ISUBC_SW = 1 ! sw: sub-col clds with prescribed seeds
!     integer, parameter :: ISUBC_SW = 2 ! sw: sub-col clds with random seeds

!     integer, parameter :: ISUBC_LW = 0 ! lw: without sub-col clds approx
!     integer, parameter :: ISUBC_LW = 1 ! lw: sub-col clds with prescribed seeds
!     integer, parameter :: ISUBC_LW = 2 ! lw: sub-col clds with random seeds
      ilons(1) = nlons
      call system_clock(count, count_rate, count_max)
      tstart = count*1.d0/count_rate
      hlw = 0; swh = 0
!$omp parallel do private(n,k,nt,i,j)
      do n=1,nlons*nlats
          ! n=i+(j-1)*nlons
          j = 1+(n-1)/nlons
          i = n-(j-1)*nlons
          do k=1,nlevs+1
             ! interface pressure, bottom to top, in cb (kPa)
             prsi(k) = pk(i,j,nlevs-k+2)/1000.
          enddo
          ! layer pressure, bottom to top, in cb (kPa)
          prsl(:) = prs(i,j,:)/1000.
          prslk = prsl
          do nt=1,ntrac
          do k=1,nlevs
             ! tracers (clipped at qmin)
             gq(k,nt) = max(qmin, tracerg(i,j,k,nt))
          enddo
          enddo
          ! sensible temp.
          gt(:) = virtempg(i,j,:)/(1.+delta*gq(:,1))
          ! omega in cb/sec
          vvel(:) = dlnpdtg(i,j,:)*prsl
!  ---  assign random seeds for sw and lw radiations
          if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
              icsdsw(1) = ixseed(i,j,1)
              icsdlw(1) = ixseed(i,j,2)
          endif
          fluxr_tmp(:) = fluxr(i,j,:)
          call grrad                                                    &
!  ---  inputs:
           ( prsi,prsl,prslk,gt,gq(1,1),gq,vvel,slmsk(i,j),             &
             lons(i,j),lats(i,j),tsea(i,j),                             &
             sheleg(i,j),sncovr(i,j),snoalb(i,j),                       &
             zorl(i,j),hprime(i,j,1),                                   &
             alvsf(i,j),alnsf(i,j),alvwf(i,j),                          &
             alnwf(i,j),facsf(i,j),facwf(i,j),                          &
             fice(i,j),tisfc(i,j),                                      &
             ! last three on next line only used for climo ozone
             solcon,coszen(i,j),coszdg(i,j),1,1,facoz,                  &
             cv(i,j),cvt(i,j),cvb(i,j),                                 &
             ! f_ice,f_rain,r_rime,flgmin only used for Ferrier microphysics
             IOVR_SW,IOVR_LW,f_ice,f_rain,r_rime,flgmin,                &
             icsdsw,icsdlw,NUM_P3D,NTCLW,NCLD,NTOZ,NTRAC,NFXR,          &
             dtlw,dtsw,lsswr,lslwr,lssav,sas_shal,norad_precip,         &
             crick_proof, ccnorm,                                       &
             1,1,nlevs,1,0,.false.,1,1,                                 &
!  ---  outputs:
             swh_tmp,sfcnsw(i,j),sfcdsw(i,j),                           &
             sfalb(i,j),                                                &
             hlw_tmp,sfcdlw(i,j),tsflw(i,j),                            &
             sfcemis(i,j),cldcov_tmp,                                   &
!  ---  input/output:
             fluxr_tmp                                                  &
             )
             ! avoid array temporaries by using these instead
             ! of passing non-contiguous slices.
             hlw(i,j,:) = hlw_tmp(:)
             swh(i,j,:) = swh_tmp(:)
             fluxr(i,j,:) = fluxr_tmp(:)
             !cldcov(i,j,:) = cldcov_tmp(:)

   enddo ! loop over horiz grid points
!$omp end parallel do 
   print *,'min/max swh',minval(swh),maxval(swh)
   print *,'min/max hlw',minval(hlw),maxval(hlw)
   call system_clock(count, count_rate, count_max)
   tend = count*1.d0/count_rate
   print *,'time in grrad = ',tend-tstart
!$omp parallel
   i = omp_get_num_threads()
!$omp end parallel
   open(7,file='grrad.dat',form='unformatted')
   if (i > 1) then
      read(7) cldcov
      print *,'sw diff',maxval(abs(cldcov-swh))
      read(7) cldcov
      print *,'lw diff',maxval(abs(cldcov-hlw))
   else
      write(7) swh
      write(7) hlw
   endif
   close(7)
   stop
   endif ! call radiation this time step?

   if (.not. newsas .or. ras) then
     print *,'old sas and ras not yet supported...'
     stop
   end if
   nnrcm=1 ! random numbers only used for old sas and ras
   ncw=1 ! only used for Ferrier microphysics (not yet supported)
   clstp=1.0 ! legacy parameter, not used
! hour of day at beginning of previous time step.
   solhr = mod(fhour - dt + idate_start(1), 24.0) 
! interpolate oz forcing to model latitudes, day of year.
   ozplout = 0.
   if (ntoz .gt. 0) then
    call ozinterpol(nlats,idate_start,fhour, &
    ozjindx1,ozjindx2,ozplin,ozplout,ozddy)
   endif
   ! random number needed for RAS and old SAS
   rann(:) = 0.6
   fscav = 0 ! only relevant for RAS

   call system_clock(count, count_rate, count_max)
   tstart = count*1.d0/count_rate
! physics loop over horiz. grid point.
!$omp parallel do private(n,k,nt,i,j,&
!$omp &gt,gu,gv,vvel,prsi,prsl,prsik,prslk,&
!$omp &phii,phil,dphi,ozp,gq,dpshc,dum1,&
!$omp &adt,adu,adv,adq,slc_tmp,stc_tmp,smc_tmp,&
!$omp &phy3d,phy2d,hlw_tmp,swh_tmp,hprime_tmp,&
!$omp &upd_mf,dwn_mf,det_mf,dkh,rnp,&
!$omp &acv,acvt,acvb,rqtk,&
!$omp &dt3dt,dq3dt,du3dt,dv3dt)
   do n=1,nlons*nlats
      ! n=i+(j-1)*nlons
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      do nt=1,ntrac
      do k=1,nlevs
         ! tracers (clipped at qmin)
         gq(k,nt) = max(qmin, tracerg(i,j,k,nt))
      enddo
      enddo
      ! sensible temp.
      gt(:) = virtempg(i,j,:)/(1.+delta*gq(:,1))
      ! omega in Pa/sec
      vvel(:) = dlnpdtg(i,j,:)*prs(i,j,:)
      ! horizontal winds.
      gu(:) = ug(i,j,:); gv(:) = vg(i,j,:)
      do k=1,nlevs+1
         ! interface pressure, bottom to top, in Pa
         prsi(k) = pk(i,j,nlevs-k+2)
      enddo
      ! layer pressure, bottom to top, in Pa.
      prsl(:) = prs(i,j,:)
      prslk = prsl**rk
      prsik = prsi**rk
      phii(1) = 0. ! topographic height not included
      do k=1,nlevs
         dphi = (prsi(k)-prsi(k+1))*(rd*virtempg(i,j,k))/(prsi(k) + prsi(k+1))
         phil(k) = phii(k) + dphi
         phii(k+1) = phil(k) + dphi
      enddo
      ozp(:,:) = ozplout(:,j,:)
      dpshc(1)     = 0.3  * prsi(1)
      phy3d(:,:) = phy_f3d(i,j,:,:)
      phy2d(:) = phy_f2d(i,j,:)
      hprime_tmp(:) = hprime(i,j,:)
      swh_tmp(:) = swh(i,j,:); hlw_tmp(:) = hlw(i,j,:)
      slc_tmp(:) = slc(i,j,:); smc_tmp(:) = smc(i,j,:); stc_tmp(:) = stc(i,j,:)
!  ====================  definition of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     ix, im   - integer, horiz dimention and num of used pts      1    !
!     levs     - integer, vertical layer dimension                 1    !
!     lsoil    - integer, number of soil layers                    1    !
!     lsm      - integer, flag for land surface model to use       1    !
!                =0  for osu lsm; =1  for noah lsm                      !
!     ntrac    - integer, number of tracers                        1    !
!     ncld     - integer, number of cloud species                  1    !
!     ntoz     - integer, ozone location in the tracer array       1    !
!     ntclw     - integer, cloud condensate location in the tracer  1    !
!                         array                                    1    !
!     nmtvr    - integer, number of topographic variables such as  1    !
!                         variance etc used in the GWD parameterization !
!     nrcm     - integer, second dimension for the random number   1    !
!                         array rann                                    !
!     ko3      - integer, number of layers for ozone data          1    !
!     lonf,latg- integer, number of lon/lat points                 1    !
!     jcap     - integer, number of spectral wave truncation       1    !
!                         used only by sascnv shalcnv                   !
!     num_p3d  - integer, number of 3D arrays needed for           1    !
!                          microphysics                                 !
!     num_p2d  - integer, number of 2D arrays needed for           1    !
!                         microphysics                                  !
!     kdt       -integer, number of the current time step          1    !
!     lat       -integer, latitude index - used for debug prints   1    !
!     me        -integer, pe number - used for debug prints        1    !
!     pl_coeff - integer, number coefficients in ozone forcing     1    !
!     nlons    - integer, number of total grid points in a latitude     !
!                         circle through a point                   im   !
!     ncw      - integer, range of droplet number concentrations for    !
!                         Ferrier microphysics                     2    !
!     flgmin   - real, range of  minimum large ice fraction for         !
!                         Ferrier microphys                        2    !
!     crtrh    - real, critical relative humidity at the surface, PBL   !
!                      top and at the top of the atmosphere        3    !
!     cdmbgwd  - real, multiplication factors for cdmb and gwd     2    !
!     ccwf     - real, multiplication factor for critical cloud         !
!                      workfunction for RAS                        2    !
!     dlqf     - real, factor for cloud condensate detrainment from     !
!                      cloud edges (RAS)                           2    !
!     ctei_rm  - real, critical cloud top entrainment instability  2    !
!                      criteria (used if mstrat=.true.)                 !
!     clstp    - real, index used by cnvc90 (for convective clouds)1    !
!                      legacy stuff - does not affect forecast          !
!     dtp      - real, physics time step in seconds                1    !
!     dtf      - real, dynamics time step in seconds               1    !
!     fhour    - real, forecast hour                               1    !
!     solhr    - real, fcst hour at the end of prev time step      1    !
!     slag     - real, equation of time ( radian )                 1    !
!     sdec,cdec- real, sin and cos of the solar declination angle  1    !
!     sinlat   - real, sin of latitude                             im   !
!     coslat   - real, cos of latitude                             im   !
!     pgr      - real, surface pressure (Pa)                       im   !
!     ugrs,vgrs- real, u/v component of layer wind              ix,levs !
!     tgrs     - real, layer mean temperature ( k )             ix,levs !
!     qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
!     vvel     - real, layer mean vertical velocity (Pa/s)      ix,levs !
!     prsi     - real, pressure at layer interfaces             ix,levs+1
!     prsl     - real, mean layer presure                       ix,levs !
!     prsik    - real, Exner function at layer interface        ix,levs+1
!     prslk    - real, Exner function at layer                  ix,levs !
!     phii     - real, interface geopotential height            ix,levs+1
!     phil     - real, layer geopotential height                ix,levs !
!     rann     - real, random number array (0-1)                ix,nrcm !
!     prdout   - real, ozone forcing data                       ix,ko3,pl_coeff!
!     poz      - real, ozone forcing data level pressure (ln(Pa))  ko3  !
!     dpshc    - real, maximum pressure depth for shallow convection im   !
!     hprime   - real, orographic std dev                       ix,nmtvr!
!     xlon,xlat- real, longitude and latitude ( radian )           im   !
!     slope    - real, sfc slope type for lsm                      im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, max snow albedo over land (for deep snow)   im   !
!     tg3      - real, deep soil temperature                       im   !
!     slmsk    - real, sea/land/ice mask (=0/1/2)                  im   !
!     vfrac    - real, vegetation fraction                         im   !
!     vtype    - real, vegetation type                             im   !
!     stype    - real, soil type                                   im   !
!     uustar   - real, boundary layer parameter                    im   !
!     oro      - real, orography                                   im   !
!     oro_uf   - real, unfiltered orography                        im   !
!     coszen   - real, avg cosz over daytime sw radiation interval im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdlw   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     tsflw    - real, sfc air (layer 1) temp over lw interval (k) im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     sfalb    - real, mean sfc diffused sw albedo                 im   !
!     swh      - real, total sky sw heating rates ( k/s )       ix,levs !
!     hlw      - real, total sky lw heating rates ( k/s )       ix,levs !
!     ras      - logical, flag for ras convection scheme           1    !
!     pre_rad  - logical, flag for testing purpose                 1    !
!     ldiag3d  - logical, flag for 3d diagnostic fields            1    !
!     lggfs3d  - logical, flag for 3d diagnostic fields for gocart 1    !
!     lssav    - logical, flag controls data store and output      1    !
!     lssav_cc - logical, flag for save data for ocean coupling    1    !
!     flipv    - logical, flag for vertical direction flip (ras)   1    !
!     xkzm_m   - real, background vertical diffusion for momentum  1    !
!     xkzm_h   - real, background vertical diffusion for heat, q   1    !
!     xkzm_s   - real, sigma threshold for background mom. diffusn 1    !
!     psautco  - real, auto conversion coeff from ice to snow      2    !
!     prautco  - real, auto conversion coeff from cloud to rain    2    !
!     evpco    - real, coeff for evaporation of largescale rain    1    !
!     wminco   - real, water and ice minimum threshold for Zhao    1    !
!     old_monin- logical, flag for diff monin schemes              1    !
!     cnvgwd   - logical, flag for conv gravity wave drag          1    !
!     shal_cnv - logical, flag for calling shallow convection      1    !
!     sashal   - logical, flag for new shallow conv scheme         1    !
!     newsas   - logical, flag for new sas conv scheme             1    !
!     cal_pre  - logical, flag controls precip type algorithm      1    !
!     mom4ice  - logical, flag controls mom4 sea-ice               1    !
!     mstrat   - logical, flag for moorthi approach for stratus    1    !
!     trans_trac-logical, flag for convective transport of tracers 1    !
!     nst_fcst  -integer, flag 0 for no nst, 1 for uncoupled nst        !
!                          and 2 for coupled NST                   1    !
!     moist_adj- logical, flag for moist convective adjustment     1    !
!     fscav    - real, tracer convective scavenging coefficient ntrac-ncld-1!
!     thermodyn_id - integer, valid for GFS only for get_prs/phi   1    !
!     sfcpress_id  - integer, valid for GFS only for get_prs/phi   1    !
!     gen_coord_hybrid - logical for Henry's gen coord             1    !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tisfc    - real, sea-ice temperature                         im   !
!     tsea     - real, ground surface temperature ( k )            im   !
!     tprcp    - real, total precipitation                         im   !
!     the following three variables do not affect the forecast          !
!     cv,       -real, convective clouds amountt                   im   !
!     cvb       -real, convective clouds base pressure (kPa)       im   !
!     cvt       -real, convective clouds top  pressure (kPa)       im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     snwdph   - real, actual snow depth (mm) over land/sea ice    im   !
!     weasd    - real, water equiv of accumulated  snow depth (kg/m**2)
!                      over land and sea ice                       im   !
!     sncovr   - real, snow cover over land                        im   !
!     zorl     - real, surface roughness                           im   !
!     canopy   - real, canopy water                                im   !
!     ffmm     - real, fm parameter from PBL scheme                im   !
!     ffhh     - real, fh parameter from PBL scheme                im   !
!     f10m     - real, fm at 10m                                   im   !
!     srunoff  - real, surface water runoff (from lsm)             im   !
!     evbsa    - real, noah lsm diagnostics                        im   !
!     evcwa    - real, noah lsm diagnostics                        im   !
!     snohfa   - real, noah lsm diagnostics                        im   !
!     transa   - real, noah lsm diagnostics                        im   !
!     sbsnoa   - real, noah lsm diagnostics                        im   !
!     snowca   - real, noah lsm diagnostics                        im   !
!     soilm    - real, soil moisture                               im   !
!     tmpmin   - real, min temperature at 2m height (k)            im   !
!     tmpmax   - real, max temperature at 2m height (k)            im   !
!     dusfc    - real, u component of surface stress               im   !
!     dvsfc    - real, v component of surface stress               im   !
!     dtsfc    - real, sensible heat flux (w/m2)                   im   !
!     dqsfc    - real, latent heat flux (w/m2)                     im   !
!     totprcp  - real, accumulated total precipitation (kg/m2)     im   !
!     gflux    - real, groud conductive heat flux                  im   !
!     dlwsfc   - real, time accumulated sfc dn lw flux ( w/m**2 )  im   !
!     ulwsfc   - real, time accumulated sfc up lw flux ( w/m**2 )  im   !
!     suntim   - real, sunshine duration time (s)                  im   !
!     runoff   - real, total water runoff                          im   !
!     ep       - real, potential evaporation                       im   !
!     cldwrk   - real, cloud workfunction (valid only with sas)    im   !
!     dugwd    - real, vertically integrated u change by OGWD      im   !
!     dvgwd    - real, vertically integrated v change by OGWD      im   !
!     psmean   - real, surface pressure (kPa)                      im   !
!     cnvprcp  - real, accumulated convective precipitation (kg/m2)im   !
!     spfhmin  - real, minimum specific humidity                   im   !
!     spfhmax  - real, maximum specific humidity                   im   !
!     dt3dt    - real, temperature change due to physics           ix,levs,6 !
!     dq3dt    - real, moisture change due to physics              ix,levs,5+pl_coeff!
!     du3dt    - real, u momentum change due to physics            ix,levs,4 !
!     dv3dt    - real, v momentum change due to physics            ix,levs,4 !
!     acv      - real,  array containing accumulated convective clouds im   !
!     acvb,acvt- real,  arrays used by cnvc90                      im   !
!     slc      - real, liquid soil moisture                     ix,lsoil!
!     smc      - real, total soil moisture                      ix,lsoil!
!     stc      - real, soil temperature                         ix,lsoil!
!     upd_mf   - real, convective updraft mass flux             ix,levs !
!     dwn_mf   - real, convective downdraft mass flux           ix,levs !
!     det_mf   - real, convective detrainment mass flux         ix,levs !
!     dkh      - real, vertical diffusion coefficient (gocart)  ix,levs !
!     rnp      - real, n cloud precipitation rate     (gocart)  ix,levs !
!     phy_f3d  - real, 3d arrays saved for restart              ix,levs,num_p3d!
!     phy_f2d  - real, 2d arrays save for restart               ix,num_p2d!
!     dlwsfc_cc- real, sfc dnwd lw flux (w/m**2) for ocn coupling  im   !
!     ulwsfc_cc- real, sfc upwd lw flux (w/m**2) for ocn coupling  im   !
!     swsfc_cc - real, sfc net sw  flux (w/m**2) for ocn coupling  im   !
!     dusfc_cc - real, sfc u-wind                for ocn coupling  im   !
!     dvsfc_cc - real, sfc v-wind                for ocn coupling  im   !
!     dtsfc_cc - real, sfc sensible heat flux    for ocn coupling  im   !
!     dqsfc_cc - real, sfc moisture flux (evap)  for ocn coupling  im   !
!     precr_cc - real, total precipitation       for ocn coupling  im   !
!
!     xt       - real, heat content in DTL                         im   !
!     xs       - real, salinity  content in DTL                    im   !
!     xu       - real, u-current content in DTL                    im   !
!     xv       - real, v-current content in DTL                    im   !
!     xz       - real, DTL thickness                               im   !
!     zm       - real, MXL thickness                               im   !
!     xtts     - real, d(xt)/d(ts)                                 im   !
!     xzts     - real, d(xz)/d(ts)                                 im   !
!     d_conv   - real, thickness of Free Convection Layer (FCL)    im   !
!     ifd      - real, index to start DTM run or not               im   !
!     dt_cool  - real, Sub-layer cooling amount                    im   !
!     Qrain    - real, sensible heat flux due to rainfall (watts)  im   !
!                                                                       !
!  outputs:                                                             !
!     adt      - real, updated temperature                        ix,levs !
!     adq      - real, updated tracers                            ix,levs,ntrac!
!     adu      - real, updated zonal wind                         ix,levs !
!     adv      - real, update meridional wind                     ix,levs !
!     t2m,q2m  - real, 2 meter temperature and humidity            im   !
!     u10m,v10m- real, 10 meater u/v wind speed                    im   !
!     zlvl     - real, layer 1 height (m)                          im   !
!     psurf    - real, surface pressure (Pa)                       im   !
!     hpbl     - real, pbl height (m)                              im   !
!     pwat     - real, precipitable water                          im   !
!     t1       - real, layer 1 temperature (K)                     im   !
!     q1       - real, layer 1 specific humidity (kg/kg)           im   !
!     u1       - real, layer 1 zonal wind (m/s)                    im   !
!     v1       - real, layer 1 merdional wind (m/s)                im   !
!     chh      - real, thermal exchange coefficient                im   !
!     cmm      - real, momentum exchange coefficient               im   !
!     dlwsfci  - real, instantaneous sfc dnwd lw flux ( w/m**2 )   im   !
!     ulwsfci  - real, instantaneous sfc upwd lw flux ( w/m**2 )   im   !
!     dswsfci  - real, instantaneous sfc dnwd sw flux ( w/m**2 )   im   !
!     uswsfci  - real, instantaneous sfc upwd sw flux ( w/m**2 )   im   !
!     dtsfci   - real, instantaneous sfc sensible heat flux        im   !
!     dqsfci   - real, instantaneous sfc latent heat flux          im   !
!     gfluxi   - real, instantaneous sfc ground heat flux          im   !
!     epi      - real, instantaneous sfc potential evaporation     im   !
!     smcwlt2  - real, wilting point (volumetric)                  im   !
!     smcref2  - real, soil moisture threshold (volumetric)        im   !
!
!     gsoil    - real                                              im   !
!     gtmp2m   - real                                              im   !
!     gustar   - real                                              im   !
!     gpblh    - real                                              im   !
!     gu10m    - real                                              im   !
!     gv10m    - real                                              im   !
!     gzorl    - real                                              im   !
!     goro     - real                                              im   !
!   
!     xmu_cc   - real, cosine of zenith angle at time step         im   !
!     dlw_cc   - real, sfc dnwd lw flux at time step for ocn cpl   im   !
!     dsw_cc   - real, sfc dnwd sw flux at time step for ocn cpl   im   !
!     snw_cc   - real, lower atms snow fall rate for ocn cpl       im   !
!     lprec_cc - real, lower atms rain fall rate for ocn cpl       im   !

!     tref     - real, Reference Temperature                       im   !
!     z_c      - real, Sub-layer cooling thickness                 im   !
!     c_0      - real, coefficient1 to calculate d(Tz)/d(Ts)       im   !
!     c_d      - real, coefficient2 to calculate d(Tz)/d(Ts)       im   !
!     w_0      - real, coefficient3 to calculate d(Tz)/d(Ts)       im   !
!     w_d      - real, coefficient4 to calculate d(Tz)/d(Ts)       im   !
!     rqtk     - real, mass change due to moisture variation       im   !
            call gbphys                                                 &
!  ---  inputs:
     &    ( 1,1,nlevs,lsoil,1,ntrac,ncld,ntoz,ntclw,            &
     &      nmtvr,nnrcm,levozp,nlons,nlats,ntrunc,num_p3d,num_p2d,          &
     &      nstep,1,0,pl_coeff,ilons,ncw,flgmin,crtrh,cdmbgwd,       &
     &      ccwf,dlqf,ctei_rm,clstp,dt,dt,fhour,solhr,                &
     &      slag,sdec,cdec,sin(lats(i,j)),cos(lats(i,j)),psg(i,j),gu,gv, &
     &      gt,gq,vvel,prsi,prsl,prslk,prsik,phii,phil,                 &
     &      rann,ozp,pl_pres,dpshc,                    &
     &      hprime_tmp,lons(i,j),lats(i,j),                      &
     &      slope (i,j),    shdmin(i,j),        &
     &      shdmax(i,j),    snoalb(i,j),        &
     &      tg3   (i,j),    slmsk (i,j),        &
     &      vfrac (i,j),    vtype (i,j),        &
     &      stype (i,j),    uustar(i,j),        &
     &      oro   (i,j),    oro_uf(i,j),        &   
     &      coszen(i,j),                                    &
     &      sfcdsw(i,j),    sfcnsw(i,j),        &
     &      sfcdlw(i,j),    tsflw (i,j),        &
     &      sfcemis(i,j),   sfalb(i,j),                 &
     &      swh_tmp,     hlw_tmp,                   &
     &      ras,pre_rad,ldiag3d,lggfs3d,lssav,lssav_cc,                 &
     &      bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,psautco,prautco, evpco, &
     &      wminco,                                                     &
     &      flipv,old_monin,cnvgwd,shal_cnv,sashal,newsas,cal_pre,      &
     &      mom4ice,mstrat,trans_trac,nst_fcst,moist_adj,fscav,         &
!    &      thermodyn_id, sfcpress_id, gen_coord_hybrid,                &
     &      0,  0, .false.,                &
!  ---  input/outputs:
     &      hice  (i,j),    fice  (i,j),        &
     &      tisfc (i,j),    tsea  (i,j),        &
     &      tprcp (i,j),    cv    (i,j),        &
     &      cvb   (i,j),    cvt   (i,j),        &
     &      srflag(i,j),    snwdph(i,j),        &
     &      sheleg(i,j),    sncovr(i,j),        &
     &      zorl  (i,j),    canopy(i,j),        &
     &      ffmm  (i,j),    ffhh  (i,j),        &
     &      f10m  (i,j),    srunoff(i,j),       &
     &      evbsa (i,j),    evcwa (i,j),        &
     &      snohfa(i,j),    transa(i,j),        &
     &      sbsnoa(i,j),    snowca(i,j),        &
     &      soilm (i,j),    tmpmin(i,j),        &
     &      tmpmax(i,j),    dusfc (i,j),        &
     &      dvsfc (i,j),    dtsfc (i,j),        &
     &      dqsfc (i,j),    geshem(i,j),        &
     &      gflux (i,j),    dlwsfc(i,j),        & 
     &      ulwsfc(i,j),    suntim(i,j),        &
     &      runoff(i,j),    ep    (i,j),        &
     &      cldwrk(i,j),    dugwd (i,j),        &
     &      dvgwd (i,j),    psmean(i,j),        &
     &      bengsh(i,j),    spfhmin(i,j),       &
     &      spfhmax(i,j),                                   &
     &      dt3dt, dq3dt, du3dt, dv3dt,                         &
     &      acv, acvb, acvt,                 &
     &      slc_tmp,smc_tmp,stc_tmp,                           &
     &      upd_mf, dwn_mf, det_mf, dkh, rnp,                 &
     &      phy3d,phy2d,                                        &
!    &      DLWSFC_cc(i,j),  ULWSFC_cc(i,j),                    &
!    &      DTSFC_cc(i,j),   SWSFC_cc(i,j),                     &
!    &      DUSFC_cc(i,j),   DVSFC_cc(i,j),                     &
!    &      DQSFC_cc(i,j),   PRECR_cc(i,j),                     &
!    &      nst_fld%xt(i,j),        nst_fld%xs(i,j),            &
!    &      nst_fld%xu(i,j),        nst_fld%xv(i,j),            &
!    &      nst_fld%xz(i,j),        nst_fld%zm(i,j),            &
!    &      nst_fld%xtts(i,j),      nst_fld%xzts(i,j),          &
!    &      nst_fld%d_conv(i,j),    nst_fld%ifd(i,j),           &
!    &      nst_fld%dt_cool(i,j),   nst_fld%Qrain(i,j),         &
     &      dum1,dum1,dum1,dum1,dum1,&
     &      dum1,dum1,dum1,dum1,dum1,&
     &      dum1,dum1,dum1,dum1,dum1,&
     &      dum1,dum1,dum1,dum1,dum1,&
!  ---  outputs:
     &      adt, adq, adu, adv,                                         &
     &      t2m   (i,j),    q2m   (i,j),        &
     &      u10m  (i,j),    v10m  (i,j),        &
     &      zlvl  (i,j),    psurf (i,j),        &
     &      hpbl  (i,j),    pwat  (i,j),        &
     &      t1    (i,j),    q1    (i,j),        &
     &      u1    (i,j),    v1    (i,j),        &
     &      chh   (i,j),    cmm   (i,j),        &
     &      dlwsfci(i,j),   ulwsfci(i,j),       &
     &      dswsfci(i,j),   uswsfci(i,j),       &
     &      dtsfci(i,j),    dqsfci(i,j),        &
     &      gfluxi(i,j),    epi   (i,j),        &
     &      smcwlt2(i,j),   smcref2(i,j),       &
     &      gsoil(i,j),     gtmp2m(i,j),        &
     &      gustar(i,j),    gpblh(i,j),         &
     &      gu10m(i,j),     gv10m(i,j),         &
     &      gzorl(i,j),     goro(i,j),          &
!    &      XMU_cc(i,j), DLW_cc(i,j), DSW_cc(i,j),          &
!    &      SNW_cc(i,j), LPREC_cc(i,j),                         &
!    &      nst_fld%Tref(i,j),       nst_fld%z_c(i,j),          &
!    &      nst_fld%c_0 (i,j),       nst_fld%c_d(i,j),          &
!    &      nst_fld%w_0 (i,j),       nst_fld%w_d(i,j),          &
!    &      rqtk)
     &      dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,rqtk)
     ! compute tendencies, 
     ! converting sensible temp back to virt temp.
     dtdt(i,j,:) = adt*(1.+delta*adq(:,1))-gt*(1.+delta*gq(:,1))
     dtracersdt(i,j,:,:) = adq-gq
     dudt(i,j,:) = adu-gu
     dvdt(i,j,:) = adv-gv
     dpsdt(i,j) = rqtk(1) ! not used
     ! avoid array temporaries by using these instead
     ! of passing non-contiguous slices.
     phy_f3d(i,j,:,:) = phy3d(:,:)
     phy_f2d(i,j,:) = phy2d(:)
     smc(i,j,:) = smc_tmp(:)
     slc(i,j,:) = slc_tmp(:)
     stc(i,j,:) = stc_tmp(:)
   enddo ! end loop over horiz grid points
!$omp end parallel do 
   call system_clock(count, count_rate, count_max)
   tend = count*1.d0/count_rate
   print *,'time in gbphys = ',tend-tstart
   print *,minval(dtdt),maxval(dtdt)
   print *,minval(dudt),maxval(dudt)
   print *,minval(dvdt),maxval(dvdt)
   print *,minval(dtracersdt(:,:,:,1)),maxval(dtracersdt(:,:,:,1))
   print *,minval(dtracersdt(:,:,:,2)),maxval(dtracersdt(:,:,:,2))
   print *,minval(dtracersdt(:,:,:,3)),maxval(dtracersdt(:,:,:,3))
   stop

! compute physics tendencies in spectral space
!$omp parallel do private(k,nt)
   do k=1,nlevs
      call grdtospec(dtdt(:,:,k), dvirtempspecdt(:,k))
      do nt=1,ntrac
         call grdtospec(dtracersdt(:,:,k,nt), dtracerspecdt(:,k,nt))
      enddo
      call getvrtdivspec(dudt(:,:,k),dvdt(:,:,k),ddivspecdt(:,k),dvrtspecdt(:,k),rerth)
   enddo
!$omp end parallel do 
   !call grdtospec(dpsdt, dlnpsspecdt)
   dlnpsspecdt=0 ! physics does not change surface pressure
   dvirtempspecdt = dvirtempspecdt/dt
   dvrtspecdt = dvrtspecdt/dt
   ddivspecdt = ddivspecdt/dt
   dtracerspecdt = dtracerspecdt/dt

   deallocate(ozplout)
   deallocate(coszdg)
   !deallocate(cldcov)
   deallocate(dtdt,dudt,dvdt)
   deallocate(dtracersdt)
   deallocate(dpsdt)

   return
 end subroutine getphytend

 subroutine getvaliddate(fhour,idate_start,id,jd)
    ! Compute valid time from initial date and forecast hour
    ! (using NCEP w3lib)
    real(r_kind), intent(in) :: fhour
    integer, dimension(8), intent(out) :: id,jd
    real(r_kind), dimension(5):: fh
    integer, intent(in),  dimension(4) :: idate_start
    fh=0; id=0; jd=0
    fh(2)=fhour    ! relative time interval in hours
    id(1)=idate_start(4) ! year
    id(2)=idate_start(2) ! month
    id(3)=idate_start(3) ! day
    id(4)=0        ! time zone
    id(5)=idate_start(1) ! hour
    call w3movdat(fh,id,jd)
    !     JDAT       INTEGER NCEP ABSOLUTE DATE AND TIME
    !                (YEAR, MONTH, DAY, TIME ZONE,
    !                 HOUR, MINUTE, SECOND, MILLISECOND)
    !idate_valid(1)=jd(5) ! hour
    !idate_valid(2)=jd(2) ! mon
    !idate_valid(3)=jd(3) ! day
    !idate_valid(4)=jd(1) ! year
    return
 end subroutine getvaliddate

  SUBROUTINE ozinterpol(nlats,IDATE,FHOUR,&
      jindx1,jindx2,ozplin,ozplout,ddy)
      implicit none
      integer  j,j1,j2,l,nc,n1,n2
      real(r_kind) fhour,tem, tx1, tx2
!
      integer  JINDX1(nlats), JINDX2(nlats)
      integer  idate(4),nlats
      integer  IDAT(8),JDAT(8)
!
      real(r_kind) ozplin(latsozp,levozp,pl_coeff,timeoz)
      real(r_kind) DDY(nlats)
      real(r_kind) ozplout(levozp,nlats,pl_coeff)
      real(r_kind) RINC(5), rjday
      integer jdow, jdoy, jday
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. PL_time(1)) RJDAY = RJDAY+365.
!
      n2 = timeoz + 1
      do j=1,timeoz
        if (rjday .lt. pl_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)     n1 = n1 + timeoz
      if (n2 > timeoz) n2 = n2 - timeoz

      tx1 = (pl_time(n2) - rjday) / (pl_time(n2) - pl_time(n1))
      tx2 = 1.0 - tx1
!
      do nc=1,pl_coeff
        DO L=1,levozp
          DO J=1,nlats
            J1  = JINDX1(J)
            J2  = JINDX2(J)
            TEM = 1.0 - DDY(J)
            ozplout(L,j,nc) = &
            tx1*(TEM*ozplin(J1,L,nc,n1)+DDY(J)*ozplin(J2,L,nc,n1)) &
          + tx2*(TEM*ozplin(J1,L,nc,n2)+DDY(J)*ozplin(J2,L,nc,n2))
          ENDDO
        ENDDO
      enddo
!
      RETURN
  end subroutine ozinterpol

end module phy_run
