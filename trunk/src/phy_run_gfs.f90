 module phy_run
! compute physics tendencies for gfs physics
! getphytend: compute total physics tendencies.
! Public subroutines:
! getphytend: compute tendencies in spectral space.

 use params, only: nlevs,nlons,nlats,ntrunc,ndimspec,ntrac,nmtvr,idate_start,dt,&
 fhzer,fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,&
 ncw,iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,&
 timestepsperhr,psautco,mstrat,pre_rad,bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,ntoz,ntclw
 use kinds, only: r_kind,r_single,r_double
 use shtns, only: grdtospec, getvrtdivspec, lons, lats
 use grid_data, only: virtempg,dlnpdtg,tracerg,ug,vg
 use pressure_data, only:  prs,psg,pk,ak,bk
 use phy_data, only: flx_init,solcon,slag,sdec,cdec,nfxr,ncld,bfilt,&
    lsoil,timeoz,latsozp,levozp,pl_coeff,ozplin,pl_pres,pl_time,&
    slmsk,sheleg,sncovr,snoalb,zorl,hprime,alvsf,ozjindx1,ozjindx2,ozddy,&
    alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,coszen,cv,cvt,cvb,sfcemis,&
    sfcnsw,sfcdsw,sfalb,sfcdlw,tsflw,slc,smc,stc,slope,shdmin,shdmax, &
    vfrac,tg3,stype,vtype,oro,oro_uf,uustar,phy_f3d,phy_f2d,&
    pl_lat,cldcov,fluxr,hlw,swh, &
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
    spfhmin,  spfhmax,        &
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
    gzorl,     goro,          &
    bengsh
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
               eps => con_eps, omega => con_omega, cvap => con_cvap, &
               grav => con_g, pi => con_pi, fv => con_fvirt, rk => con_rocp

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
   real(r_double), intent(in) :: t,dtx
   real(r_kind) dtp
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   real(r_kind), parameter :: qmin = 1.e-10 ! min value for clipping tracers
   real(r_kind), parameter :: typical_pgr = 95000.0
   real(r_kind)  :: fscav(ntrac-ncld-1),&
   fhour,dtsw,dtlw,facoz,clstp,solhr,dphi,dpshc(1),&
   gt(nlevs),prsl(nlevs),prsi(nlevs+1),vvel(nlevs),&
   f_ice(nlevs),f_rain(nlevs),r_rime(nlevs),&
   prslk(nlevs),gq(nlevs,ntrac),prsik(nlevs+1),si_loc(nlevs+1),phii(nlevs+1),&
   phil(nlevs),rann(1),acv(1),acvb(1),acvt(1),adt(nlevs),adu(nlevs),adv(nlevs),&
   dum1(1),rqtk(1),upd_mf(nlevs),dwn_mf(nlevs),det_mf(nlevs),dkh(nlevs),rnp(nlevs),&
   adq(nlevs,ntrac),dt3dt(nlevs,6),du3dt(nlevs,4),dv3dt(nlevs,4),dq3dt(nlevs,5+pl_coeff),&
   phy3d(nlevs,num_p3d),phy2d(num_p2d),hlw_tmp(nlevs),swh_tmp(nlevs),fluxr_tmp(nfxr),&
   cldcov_tmp(nlevs),hprime_tmp(nmtvr),slc_tmp(lsoil),smc_tmp(lsoil),stc_tmp(lsoil),&
   gu(nlevs),gv(nlevs)
   integer :: icsdsw(1),icsdlw(1),ilons(1), ipsd0
   real(r_kind), allocatable, dimension(:,:) :: coszdg,dpsdt
   real(r_kind), allocatable, dimension(:,:,:) :: ozplout,dtdt,dudt,dvdt
   real(4), allocatable, dimension(:,:,:) :: work4
   real(r_kind), allocatable, dimension(:,:,:,:) :: dtracersdt
   logical :: sas_shal ! sas shallow convection.
   ! set to true on first call to indicate
   ! that random number seed for mcica radiation needs to be initialized.
   logical, save :: first=.true. 
   logical :: lsswr=.false. ! sw rad call
   logical :: lslwr=.false. ! lw rad call
   logical :: lszer=.false. ! zero out flux accum.
   ! flag for saving physics data in gbphys
   ! (always true, except in digitial filter)
   logical :: lssav=.true. 
   logical :: lssav_cc=.false. ! flag for coupling
   logical :: flipv=.true. ! vert dir flip for RAS
   logical :: ldiag3d = .false.
   logical :: lggfs3d = .false.
   integer :: idat(8), jdat(8)
   type (random_stat) :: stat
   integer :: numrdm(nlons*nlats*2), ixseed(nlons,nlats,2)
   integer :: i,j,k,n,nt,ipseed,nstep,nswr,nlwr,nnrcm,nszer
   integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds
   real(8) tstart,tend
   integer(8) count, count_rate, count_max
   logical :: testomp=.true.  ! openmp debug flag

   dtp = dtx

   allocate(ozplout(levozp,pl_coeff,nlats))
   allocate(coszdg(nlons,nlats))
   allocate(dtdt(nlons,nlats,nlevs),dudt(nlons,nlats,nlevs),dvdt(nlons,nlats,nlevs))
   allocate(dtracersdt(nlons,nlats,nlevs,ntrac))
   allocate(dpsdt(nlons,nlats))

! is it a radiation time step (long or short wave)?
   fhour = t/3600.
   nswr = int(fhswr*timestepsperhr)
   nlwr = int(fhlwr*timestepsperhr)
   nszer = int(fhzer*timestepsperhr)
   if (t > 0) then
      nstep = nint(t/dt)
      if (nswr .eq. 0) then
         lsswr = .true.
      else
         lsswr = mod(nstep,nswr) .eq. 0
      endif
      if (nlwr .eq. 0) then
         lslwr = .true.
      else
         lslwr = mod(nstep,nlwr) .eq. 0
      endif
      lszer = mod(nstep,nszer) .eq. 0
   else
      nstep=0
      lsswr=.true.
      lslwr=.true.
      lszer=.false.
   endif
! compute forecast valid time (jdat).
   call getvaliddate(fhour,idate_start,idat,jdat)

   sas_shal = sashal .and. (.not. ras)
! get a sigma distribution for radiation-cloud initialization
   si_loc(nlevs+1)= ak(1)/typical_pgr+bk(1)
   do k=1,nlevs
      si_loc(nlevs+1-k)= ak(k+1)/typical_pgr + bk(k+1)
   enddo
! initialize on first call
   if (first) then
! generate initial permutation seed for random number generator
     if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
        ipsd0 = 17*idate_start(1) + 43*idate_start(2) + &
                37*idate_start(3) + 23*idate_start(4)
        print *,'  Radiation sub-cloud initial seed =',ipsd0,&
                ' idate =',idate_start
     endif
     first = .false.
   endif         ! end_if_first

   ! reset accumulated arrays
   if (lszer) then
      fluxr = 0.
      print *,'zeroing fluxes fh,nt=',fhour,nt,nszer
      call flx_init()
   endif
!
!===> *** ...  radiation initialization
!
   if (lsswr .or. lslwr) then ! call radiation

! the following block of code is adapted from gloopr.f
   if (fhswr .lt. dtx) then
      dtsw = dtx
   else
      dtsw  = 3600.0 * fhswr
   endif
   if (fhlwr .lt. dtx) then
      dtlw = dtx
   else
      dtlw  = 3600.0 * fhlwr
   endif
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
       ( nlons, nlats, lons, lats,                                   &
         fhswr, jdat, lsswr,                                         &
!  ---  outputs:
         solcon, slag, sdec, cdec, coszen, coszdg                    &
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
      endif ! end random number seeds generation

   call system_clock(count, count_rate, count_max)
   tstart = count*1.d0/count_rate
!$omp parallel do private(n,k,nt,i,j,prslk,f_ice,f_rain,&
!$omp& r_rime,prsi,prsl,gt,gq,vvel,icsdsw,icsdlw,fluxr_tmp,&
!$omp& swh_tmp,hlw_tmp,cldcov_tmp) schedule(dynamic)
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
       do k=1,nlevs
          ! tracers (humidity, ozone, cloud condensate)
          ! clipped to qmin
          do nt=1,ntrac
             gq(k,nt) = max(qmin,tracerg(i,j,k,nt))
          enddo
          ! compute sensible temp.
          gt(k) = virtempg(i,j,k)/(1.+fv*gq(k,1))
       enddo
       ! omega in cb/sec
       vvel(:) = dlnpdtg(i,j,:)*prsl
!  ---  assign random seeds for sw and lw radiations
       if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
           icsdsw(1) = ixseed(i,j,1)
           icsdlw(1) = ixseed(i,j,2)
       endif
       fluxr_tmp(:) = fluxr(i,j,:)
! call GFS radiation driver
       call grrad                                                   &
!  ---  inputs:
       ( prsi,prsl,prslk,gt,gq(:,1),gq,vvel,slmsk(i,j),             &
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
        fluxr_tmp                                                   &
       )
       ! avoid array temporaries by using these instead
       ! of passing non-contiguous slices.
       hlw(i,j,:) = hlw_tmp(:)
       swh(i,j,:) = swh_tmp(:)
       fluxr(i,j,:) = fluxr_tmp(:)
       ! diagnosed 3d cloud fraction.
       cldcov(i,j,:) = cldcov_tmp(:)

   enddo ! loop over horiz grid points
!$omp end parallel do 
   if (testomp) then
      print *,'min/max swh',minval(swh),maxval(swh)
      print *,'min/max hlw',minval(hlw),maxval(hlw)
      call system_clock(count, count_rate, count_max)
      tend = count*1.d0/count_rate
      print *,'time in grrad = ',tend-tstart
!$omp parallel
      i = omp_get_num_threads()
!$omp end parallel
      allocate(work4(nlons,nlats,nlevs))
      open(7,file='grrad.dat',form='unformatted')
      if (i > 1) then
         read(7) work4
         print *,'sw diff',maxval(abs(work4-swh))
         read(7) work4
         print *,'lw diff',maxval(abs(work4-hlw))
      else
         work4 = swh
         write(7) work4
         work4 = hlw
         write(7) work4
      endif
      close(7)
   endif

   endif ! call radiation this time step?

   if (.not. newsas .or. ras) then
     print *,'old sas and ras not yet supported...'
     stop
   end if
   nnrcm=1 ! random numbers only used for old sas and ras
   ncw=1 ! only used for Ferrier microphysics (not yet supported)
   clstp=1.0 ! legacy parameter, not used
! hour of day at end of previous time step.
   solhr = mod(fhour + idate_start(1), 24.0) 
! interpolate oz forcing to model latitudes, day of year.
   ozplout = 0.
   if (ntoz .gt. 0) then
      call ozinterpol(jdat,ozplout)
   endif
   ! random number needed for RAS and old SAS
   rann(:) = 0.6
   fscav = 0 ! only relevant for RAS
   ilons(1) = nlons

   call system_clock(count, count_rate, count_max)
   tstart = count*1.d0/count_rate
! physics loop over horiz. grid points.
!$omp parallel do private(n,k,i,j,&
!$omp& gt,gu,gv,gq,vvel,prsi,prsl,prsik,prslk,&
!$omp& phii,phil,dphi,dpshc,dum1,&
!$omp& adt,adu,adv,adq,slc_tmp,stc_tmp,smc_tmp,&
!$omp& phy3d,phy2d,hlw_tmp,swh_tmp,hprime_tmp,&
!$omp& upd_mf,dwn_mf,det_mf,dkh,rnp,&
!$omp& acv,acvt,acvb,rqtk,&
!$omp& dt3dt,dq3dt,du3dt,dv3dt) schedule(dynamic)
   do n=1,nlons*nlats
      ! n=i+(j-1)*nlons
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      ! tracers (humidity, ozone, cloud condensate - don't clip)
      gq = tracerg(i,j,:,:)
      ! compute sensible temp (clip humidity in computation).
      do k=1,nlevs
         gt(k) = virtempg(i,j,k)/(1.+fv*max(qmin,gq(k,1)))
      enddo
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
      phil(nlevs) = 0. ! will be recomputed in gbphys
      !phii(1) = 0. ! topographic height not included
      !do k=1,nlevs
      !   dphi = (prsi(k)-prsi(k+1))*(rd*virtempg(i,j,k))/(prsi(k) + prsi(k+1))
      !   phil(k) = phii(k) + dphi
      !   phii(k+1) = phil(k) + dphi
      !enddo
      dpshc(1)     = 0.3  * prsi(1)
      phy3d(:,:) = phy_f3d(i,j,:,:)
      phy2d(:) = phy_f2d(i,j,:)
      hprime_tmp(:) = hprime(i,j,:)
      swh_tmp(:) = swh(i,j,:); hlw_tmp(:) = hlw(i,j,:)
      slc_tmp(:) = slc(i,j,:); smc_tmp(:) = smc(i,j,:); stc_tmp(:) = stc(i,j,:)
! call GFS physics driver
      call gbphys                                                       &
!  ---  inputs:
          ( 1,1,nlevs,lsoil,1,ntrac,ncld,ntoz,ntclw,                    &
            nmtvr,nnrcm,levozp,nlons,nlats,ntrunc,num_p3d,num_p2d,      &
            nstep,1,0,pl_coeff,ilons,ncw,flgmin,crtrh,cdmbgwd,          &
            ccwf,dlqf,ctei_rm,clstp,dtp,dtp,fhour,solhr,                &
            slag,sdec,cdec,sin(lats(i,j)),cos(lats(i,j)),psg(i,j),gu,gv,&
            gt,gq,vvel,prsi,prsl,prslk,prsik,phii,phil,                 &
            rann,ozplout(:,:,j),pl_pres,dpshc,                          &
            hprime_tmp,lons(i,j),lats(i,j),                             &
            slope (i,j),    shdmin(i,j),                                &
            shdmax(i,j),    snoalb(i,j),                                &
            tg3   (i,j),    slmsk (i,j),                                &
            vfrac (i,j),    vtype (i,j),                                &
            stype (i,j),    uustar(i,j),                                &
            oro   (i,j),    oro_uf(i,j),                                &   
            coszen(i,j),                                                &
            sfcdsw(i,j),    sfcnsw(i,j),                                &
            sfcdlw(i,j),    tsflw (i,j),                                &
            sfcemis(i,j),   sfalb(i,j),                                 &
            swh_tmp,     hlw_tmp,                                       &
            ras,pre_rad,ldiag3d,lggfs3d,lssav,lssav_cc,                 &
            bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,psautco,prautco, evpco, &
            wminco,                                                     &
            flipv,old_monin,cnvgwd,shal_cnv,sashal,newsas,cal_pre,      &
            mom4ice,mstrat,trans_trac,nst_fcst,moist_adj,fscav,         &
!           thermodyn_id, sfcpress_id, gen_coord_hybrid,                &
            0,  0, .false.,                                             &
!  ---  input/outputs:
            hice  (i,j),    fice  (i,j),                                &
            tisfc (i,j),    tsea  (i,j),                                &
            tprcp (i,j),    cv    (i,j),                                &
            cvb   (i,j),    cvt   (i,j),                                &
            srflag(i,j),    snwdph(i,j),                                &
            sheleg(i,j),    sncovr(i,j),                                &
            zorl  (i,j),    canopy(i,j),                                &
            ffmm  (i,j),    ffhh  (i,j),                                &
            f10m  (i,j),    srunoff(i,j),                               &
            evbsa (i,j),    evcwa (i,j),                                &
            snohfa(i,j),    transa(i,j),                                &
            sbsnoa(i,j),    snowca(i,j),                                &
            soilm (i,j),    tmpmin(i,j),                                &
            tmpmax(i,j),    dusfc (i,j),                                &
            dvsfc (i,j),    dtsfc (i,j),                                &
            dqsfc (i,j),    geshem(i,j),                                &
            gflux (i,j),    dlwsfc(i,j),                                & 
            ulwsfc(i,j),    suntim(i,j),                                &
            runoff(i,j),    ep    (i,j),                                &
            cldwrk(i,j),    dugwd (i,j),                                &
            dvgwd (i,j),    psmean(i,j),                                &
            bengsh(i,j),    spfhmin(i,j),                               &
            spfhmax(i,j),                                               &
            dt3dt, dq3dt, du3dt, dv3dt,                                 &
            acv, acvb, acvt,                                            &
            slc_tmp,smc_tmp,stc_tmp,                                    &
            upd_mf, dwn_mf, det_mf, dkh, rnp,                           &
            phy3d,phy2d,                                                &
! nst and coupled ocean components not implemented.
!           DLWSFC_cc(i,j),  ULWSFC_cc(i,j),                            &
!           DTSFC_cc(i,j),   SWSFC_cc(i,j),                             &
!           DUSFC_cc(i,j),   DVSFC_cc(i,j),                             &
!           DQSFC_cc(i,j),   PRECR_cc(i,j),                             &
!           nst_fld%xt(i,j),        nst_fld%xs(i,j),                    &
!           nst_fld%xu(i,j),        nst_fld%xv(i,j),                    &
!           nst_fld%xz(i,j),        nst_fld%zm(i,j),                    &
!           nst_fld%xtts(i,j),      nst_fld%xzts(i,j),                  &
!           nst_fld%d_conv(i,j),    nst_fld%ifd(i,j),                   &
!           nst_fld%dt_cool(i,j),   nst_fld%Qrain(i,j),                 &
            dum1,dum1,dum1,dum1,dum1,                                   &
            dum1,dum1,dum1,dum1,dum1,                                   &
            dum1,dum1,dum1,dum1,dum1,                                   &
            dum1,dum1,dum1,dum1,dum1,                                   &
!  --   outputs:
            adt, adq, adu, adv,                                         &
            t2m   (i,j),    q2m   (i,j),                                &
            u10m  (i,j),    v10m  (i,j),                                &
            zlvl  (i,j),    psurf (i,j),                                &
            hpbl  (i,j),    pwat  (i,j),                                &
            t1    (i,j),    q1    (i,j),                                &
            u1    (i,j),    v1    (i,j),                                &
            chh   (i,j),    cmm   (i,j),                                &
            dlwsfci(i,j),   ulwsfci(i,j),                               &
            dswsfci(i,j),   uswsfci(i,j),                               &
            dtsfci(i,j),    dqsfci(i,j),                                &
            gfluxi(i,j),    epi   (i,j),                                &
            smcwlt2(i,j),   smcref2(i,j),                               &
            gsoil(i,j),     gtmp2m(i,j),                                &
            gustar(i,j),    gpblh(i,j),                                 &
            gu10m(i,j),     gv10m(i,j),                                 &
            gzorl(i,j),     goro(i,j),                                  &
!           XMU_cc(i,j), DLW_cc(i,j), DSW_cc(i,j),                      &
!           SNW_cc(i,j), LPREC_cc(i,j),                                 &
!           nst_fld%Tref(i,j),       nst_fld%z_c(i,j),                  &
!           nst_fld%c_0 (i,j),       nst_fld%c_d(i,j),                  &
!           nst_fld%w_0 (i,j),       nst_fld%w_d(i,j),                  &
!           rqtk)
            dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,dum1,rqtk)
     ! convert sensible temp back to virt temp.
     ! (clip humidity in conversion)
     do k=1,nlevs
        adt(k) = adt(k)*(1.+fv*max(qmin,adq(k,1)))
     enddo
     ! compute tendencies, 
     ! update grid data.
     dtdt(i,j,:) = adt-virtempg(i,j,:)
     dtracersdt(i,j,:,:) = adq-tracerg(i,j,:,:)
     dudt(i,j,:) = adu-ug(i,j,:)
     dvdt(i,j,:) = adv-vg(i,j,:)
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
   !print *,'min/max dtdt',minval(dtdt),maxval(dtdt)
   !print *,'min/max dudt',minval(dudt),maxval(dudt)
   !print *,'min/max dvdt',minval(dvdt),maxval(dvdt)
   !print *,'min/max dtracer1dt',minval(dtracersdt(:,:,:,1)),maxval(dtracersdt(:,:,:,1))
   !print *,'min/max dtracer2dt',minval(dtracersdt(:,:,:,2)),maxval(dtracersdt(:,:,:,2))
   !print *,'min/max dtracer3dt',minval(dtracersdt(:,:,:,3)),maxval(dtracersdt(:,:,:,3))
   !print *,'min/max tracer1',minval(tracerg(:,:,:,1)),maxval(tracerg(:,:,:,1))
   !print *,'min/max tracer2',minval(tracerg(:,:,:,2)),maxval(tracerg(:,:,:,2))
   !print *,'min/max tracer3',minval(tracerg(:,:,:,3)),maxval(tracerg(:,:,:,3))
   if (testomp) then
      call system_clock(count, count_rate, count_max)
      tend = count*1.d0/count_rate
      print *,'time in gbphys = ',tend-tstart
      print *,'min/max dtdt',minval(dtdt),maxval(dtdt)
      print *,'min/max dudt',minval(dudt),maxval(dudt)
      print *,'min/max dvdt',minval(dvdt),maxval(dvdt)
      print *,'min/max dtracer1dt',minval(dtracersdt(:,:,:,1)),maxval(dtracersdt(:,:,:,1))
      print *,'min/max dtracer2dt',minval(dtracersdt(:,:,:,2)),maxval(dtracersdt(:,:,:,2))
      print *,'min/max dtracer3dt',minval(dtracersdt(:,:,:,3)),maxval(dtracersdt(:,:,:,3))
      print *,'min/max tracer1',minval(tracerg(:,:,:,1)),maxval(tracerg(:,:,:,1))
      print *,'min/max tracer2',minval(tracerg(:,:,:,2)),maxval(tracerg(:,:,:,2))
      print *,'min/max tracer3',minval(tracerg(:,:,:,3)),maxval(tracerg(:,:,:,3))
!$omp parallel
      i = omp_get_num_threads()
!$omp end parallel
      open(7,file='gbphys.dat',form='unformatted')
      if (i > 1) then
         read(7) work4
         print *,'dtdt diff',maxval(abs(work4-dtdt))
         read(7) work4
         print *,'dudt diff',maxval(abs(work4-dudt))
         read(7) work4
         print *,'dvdt diff',maxval(abs(work4-dvdt))
         do nt=1,ntrac
            read(7) work4
            !print *,minval(work4),maxval(work4)
            print *,nt,'dtracerdt diff',maxval(abs(work4-dtracersdt(:,:,:,nt)))
         enddo
      else
         work4 = dtdt
         write(7) work4
         work4 = dudt
         write(7) work4
         work4 = dvdt
         write(7) work4
         do nt=1,ntrac
            work4 = dtracersdt(:,:,:,nt)
            write(7) work4
         enddo
      endif
      close(7)
      deallocate(work4)
      stop
   endif

! compute physics tendencies in spectral space
!$omp parallel do private(k,nt)
   do k=1,nlevs
      call grdtospec(dtdt(:,:,k), dvirtempspecdt(:,k))
      call getvrtdivspec(dudt(:,:,k),dvdt(:,:,k),dvrtspecdt(:,k),ddivspecdt(:,k),rerth)
      do nt=1,ntrac
         call grdtospec(dtracersdt(:,:,k,nt), dtracerspecdt(:,k,nt))
      enddo
      ! apply "gloopb filter"
      dvirtempspecdt(:,k) = bfilt(:)*dvirtempspecdt(:,k)
      dvrtspecdt(:,k) = bfilt(:)*dvrtspecdt(:,k)
      ddivspecdt(:,k) = bfilt(:)*ddivspecdt(:,k)
      do nt=1,ntrac
         dtracerspecdt(:,k,nt) = bfilt(:)*dtracerspecdt(:,k,nt)
      enddo
   enddo
!$omp end parallel do 
   !call grdtospec(dpsdt, dlnpsspecdt)
   dlnpsspecdt=0 ! physics does not change surface pressure
   dvirtempspecdt = dvirtempspecdt/dtx
   dvrtspecdt = dvrtspecdt/dtx
   ddivspecdt = ddivspecdt/dtx
   dtracerspecdt = dtracerspecdt/dtx

   deallocate(ozplout)
   deallocate(coszdg)
   deallocate(dtdt,dudt,dvdt)
   deallocate(dtracersdt)
   deallocate(dpsdt)

   return
 end subroutine getphytend

 subroutine getvaliddate(fhour,idate_start,id,jd)
    ! Compute valid time from initial date and forecast hour
    ! (using NCEP w3lib)
    !     JD       INTEGER NCEP ABSOLUTE DATE AND TIME
    !              (YEAR, MONTH, DAY, TIME ZONE,
    !              HOUR, MINUTE, SECOND, MILLISECOND)
    real(r_kind), intent(in) :: fhour
    integer, dimension(8), intent(out) :: id,jd
    real(4), dimension(5):: fh
    integer, intent(in),  dimension(4) :: idate_start
    fh=0; id=0; jd=0
    fh(2)=fhour          ! relative time interval in hours
    id(1)=idate_start(4) ! year
    id(2)=idate_start(2) ! month
    id(3)=idate_start(3) ! day
    id(4)=0              ! time zone
    id(5)=idate_start(1) ! hour
    call w3movdat(fh,id,jd)
 end subroutine getvaliddate

 SUBROUTINE ozinterpol(jdat,ozplout)
   ! interpolate climatological ozone forcing in space and time.
   implicit none
   integer  j,j1,j2,l,nc,n1,n2,jdow,jdoy,jday
   real(r_kind) rjday, tx1, tx2
   integer,intent(in) :: JDAT(8)
   real(r_kind), intent(out) ::  ozplout(levozp,pl_coeff,nlats)
   jdow = 0
   jdoy = 0
   jday = 0
   call w3doxdat(jdat,jdow,jdoy,jday)
   rjday = jdoy + jdat(5) / 24.
   IF (RJDAY .LT. PL_time(1)) RJDAY = RJDAY+365.
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
   do nc=1,pl_coeff
     DO L=1,levozp
       DO J=1,nlats
         J1  = OZJINDX1(J)
         J2  = OZJINDX2(J)
         ozplout(L,nc,j) = &
         tx1*((1.-OZDDY(J))*ozplin(J1,L,nc,n1)+OZDDY(J)*ozplin(J2,L,nc,n1)) +&
         tx2*((1.-OZDDY(J))*ozplin(J1,L,nc,n2)+OZDDY(J)*ozplin(J2,L,nc,n2))
       ENDDO
     ENDDO
   enddo
 end subroutine ozinterpol

end module phy_run
