module dyn_init
! initialize dynamics (including reading in initial conditions).
! also does IO of spectral data.
! public subroutines:
! init_dyn:  allocate spectral arrays (init_specdata),
!  read in initial conditions, initalize pressure variables (init_pressdata,
!  get ak,bk from IC file header), compute gradient of orography,
!  set up spectral hyperdiffusion operator (disspec), initialize arrays
!  for semi-implicit time stepping.
! wrtout: write out spectral data.
 use kinds, only: r_kind, r_single
 use sigio_module, only: sigio_sclose,sigio_swohdc,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_sropen,sigio_srdata,sigio_axdata
 use params, only: &
 nlons,nlats,nlevs,ndimspec,ntrunc,initfile,sighead,dry,ndiss,efold,jablowill,&
 heldsuarez,explicit
 use shtns, only: shtns_init, spectogrd, grdtospec, getgrad, getvrtdivspec, lap, lats, lons
 use spectral_data, only: vrtspec,divspec,virtempspec,spfhumspec,topospec,lnpsspec,&
                          disspec,init_specdata
 use pressure_data, only: ak,bk,ck,dbk,psg,prs,init_pressdata,calc_pressdata
 use grid_data, only: init_griddata, dphisdx, dphisdy, phis, ug, vg, virtempg
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
                     omega => con_omega, grav => con_g, pi => con_pi
 use semimp_data, only: init_semimpdata

 implicit none
 private
 public :: init_dyn, wrtout

 contains

 subroutine init_dyn()
    type(sigio_data) sigdata
    real(r_kind), dimension(:,:), allocatable :: psg
    integer k,lu,iret
    ! allocate arrays
    call init_specdata()
    ! initialize spherical harmonic lib
    call shtns_init(nlons,nlats,ntrunc)
    ! read initial conditions
    lu = 7
    call sigio_srohdc(lu,trim(initfile),sighead,sigdata,iret)
    if (iret .ne. 0) then
      print *,'error reading ',trim(initfile),iret
      stop
    endif
    ! convert spectral arrays to double precision complex,
    ! re-normalize coefficients.
    call copyspecin(sigdata%ps, lnpsspec)
    call copyspecin(sigdata%hs, topospec)
    do k=1,nlevs
       call copyspecin(sigdata%z(:,k),vrtspec(:,k))
       call copyspecin(sigdata%d(:,k),divspec(:,k))
       call copyspecin(sigdata%t(:,k),virtempspec(:,k))
       call copyspecin(sigdata%q(:,k,1),spfhumspec(:,k))
    enddo
    if (dry) spfhumspec = 0.
    allocate(psg(nlons,nlats))
    ! convert to ln(ps) in Pa.
    call spectogrd(lnpsspec, psg)
    psg = 1000.*exp(psg) ! convert to Pa
    call grdtospec(log(psg), lnpsspec) ! back to spectral.
    call sigio_axdata(sigdata,iret)
    call sigio_sclose(lu,iret)
    ! initialize pressure arrays.
    call init_pressdata()
    if (sighead%idvc == 2) then ! hybrid coordinate
       do k=1,nlevs+1
          ak(k) = sighead%vcoord(nlevs+2-k,1)
          bk(k) = sighead%vcoord(nlevs+2-k,2)
       enddo
       do k=1,nlevs
          dbk(k) = bk(k+1)-bk(k)
          ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
       enddo
    else
       print *,'unknown vertical coordinate type',sighead%idvc
       stop
    end if
    call init_griddata() ! allocate gridded fields
    ! jablonowoski and williamson test case or
    ! held-suarez test case
    ! (over-ride initial conditions read from file).
    if (heldsuarez) then
       call jablowill_ics()
       spfhumspec = 0.
    endif
    if (heldsuarez) then
       call heldsuarez_ics()
       spfhumspec = 0.
    endif
    ! print out max/min values of phis,ps
    call spectogrd(grav*topospec, phis)
    print *,'min/max surface geopotential',minval(phis),maxval(phis)
    call spectogrd(lnpsspec, psg)
    psg = exp(psg) 
    print *,'min/max sfc pressure (hPa)',maxval(psg/100.),minval(psg/100.)
    deallocate(psg)
    ! compute gradient of surface orography
    call getgrad(grav*topospec, dphisdx, dphisdy, rerth)
    ! hyper-diffusion operator
    disspec(:) = -(1./efold)*(lap(:)/minval(lap))**(ndiss/2)
    ! initialize arrays for semi-implicit adjustments.
    if (.not. explicit) call init_semimpdata()
 end subroutine init_dyn

 subroutine copyspecin(rspecdata,cspecdata)
    real(r_single), intent(in) :: rspecdata(2*ndimspec)
    complex(r_kind), intent(out) :: cspecdata(ndimspec)
    integer n,nn
    nn = 1
    ! factor of sqrt(2.*pi) accounts for difference in normalization
    ! between ncep libs and shtns (which uses orthonormalized norm)
    do n=1,ndimspec
       cspecdata(n) = sqrt(2.*pi)*cmplx(rspecdata(nn),rspecdata(nn+1))
       nn = nn + 2
    enddo
 end subroutine copyspecin

 subroutine copyspecout(cspecdata, rspecdata)
    real(r_single), intent(out) :: rspecdata(2*ndimspec)
    complex(r_kind), intent(in) :: cspecdata(ndimspec)
    integer n,nn
    nn = 1
    ! factor of sqrt(2.*pi) accounts for difference in normalization
    ! between ncep libs and shtns (which uses orthonormalized norm)
    do n=1,ndimspec
       rspecdata(nn) = real(cspecdata(n))/sqrt(2.*pi)
       rspecdata(nn+1) = imag(cspecdata(n))/sqrt(2.*pi)
       nn = nn + 2
    enddo
 end subroutine copyspecout

 subroutine wrtout(fh,filename)
    real(r_kind), dimension(nlons,nlats) :: psg 
    complex(r_kind), dimension(ndimspec) :: lnpsspec_tmp
    type(sigio_data) sigdata
    real(r_kind), intent(in) :: fh
    character(len=500), intent(in) :: filename
    integer k,lu,iret,lu2
    ! convert to ln(ps) in Pa.
    call spectogrd(lnpsspec, psg)
    psg = exp(psg)/1000. ! convert to cb
    call grdtospec(log(psg), lnpsspec_tmp) ! back to spectral.
    print *,'fhour =',int(fh),', min/max psg = ',maxval(10.*psg),minval(10.*psg)
    lu = 7; lu = 8
    call sigio_srohdc(lu,trim(initfile),sighead,sigdata,iret)
    if (iret .ne. 0) then
      print *,'error reading ',trim(initfile),iret
      stop
    endif
    sighead%fhour = fh
    call copyspecout(lnpsspec_tmp, sigdata%ps)
    call copyspecout(topospec, sigdata%hs)
    do k=1,nlevs
       call copyspecout(vrtspec(:,k), sigdata%z(:,k))
       call copyspecout(divspec(:,k), sigdata%d(:,k))
       call copyspecout(virtempspec(:,k), sigdata%t(:,k))
       call copyspecout(spfhumspec(:,k), sigdata%q(:,k,1))
    enddo
    call sigio_swohdc(lu2,filename,sighead,sigdata,iret)
    if (iret .ne. 0) then
      print *,'error writing ',trim(filename),iret
      stop
    endif
    call sigio_axdata(sigdata,iret)
    call sigio_sclose(lu,iret)
    call sigio_sclose(lu2,iret)
 end subroutine wrtout

 subroutine jablowill_ics()
   ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
   ! idealized baroclinic instability test case initial conditions.
   integer k
   real(r_kind) u0,etat,eta0,t0,gamma,deltat,lonc,latc,&
        up,pertrad
   real(r_kind), dimension(nlons,nlats) :: r,x,eta,etav,zs,lnpsg
   print *,'replacing initial conds with jablonowsky and williamson test case..'
   ! zonal jet.
   psg = 1.e5
   lnpsg = log(psg)
   call grdtospec(lnpsg,lnpsspec)
   call calc_pressdata(lnpsg)
   eta0 = 0.252
   u0 = 35.
   t0 = 288.
   gamma = 0.005
   deltat = 4.8e5
   etat = 0.2
   lonc = pi/9.
   latc = 2.*pi/9.
   up = 1.
   pertrad = rerth/10.
   do k=1,nlevs
      eta = prs(:,:,nlevs-k+1)/psg
      etav = 0.5*pi*(eta-eta0)
      ug(:,:,k) = u0*cos(etav)**1.5*sin(2.*lats)**2 ! eqn 2
      ! eqns 4 and 5
      virtempg(:,:,k) = t0*eta**(rd*gamma/grav)
      where (eta < etat) 
         virtempg(:,:,k) = virtempg(:,:,k) +&
         deltat*(etat-eta)**5
      end where
      ! eqn 6
      virtempg(:,:,k) = virtempg(:,:,k) +&
      0.75*(eta*pi*u0/rd)*sin(etav)*sqrt(cos(etav))*&
      ((-2.*sin(lats)**6*(cos(lats)**2+(1./3.)) + (10./63))*2.*u0*cos(etav)**1.5+&
      ((8./5.)*cos(lats)**3*(sin(lats)**2+(2./3.)) - 0.25*pi)*rerth*omega)
      !print *,k,minval(ug(:,:,k)),maxval(ug(:,:,k)),minval(virtempg(:,:,k)),&
      !maxval(virtempg(:,:,k))
   enddo
   ! surface height (to balance wind field so surface pressure can be constant).
   eta = 1; etav = 0.5*pi*(eta-eta0)
   zs = u0*cos(etav)**1.5*((-2.*sin(lats)**6*(cos(lats)**2+(1./3.)) + (10./63))*u0*cos(etav)**1.5+&
      ((8./5.)*cos(lats)**3*(sin(lats)**2+(2./3.)) - 0.25*pi)*rerth*omega)/grav
   call grdtospec(zs, topospec) 
   !print *,minval(grav*zs),maxval(grav*zs)
   !do k=1,nlevs
   !   eta = prs(:,:,nlevs-k+1)/psg
   !   print *,k,eta(1,nlats/4),virtempg(1,nlats/4,k)
   !enddo
   !stop
   ! add perturbation
   x = sin(latc)*sin(lats) + cos(latc)*cos(lats)*cos(lons-lonc)
   r = rerth*acos(x)
   do k=1,nlevs
      ! add a zonal wind perturbation
      ug(:,:,k) = ug(:,:,k) + up*exp(-(r/pertrad)**2)
      vg(:,:,k) = 0.
      call getvrtdivspec(ug(:,:,k),vg(:,:,k),vrtspec(:,k),divspec(:,k),rerth)
      call grdtospec(virtempg(:,:,k),virtempspec(:,k))
   enddo
 end subroutine jablowill_ics

 subroutine heldsuarez_ics()
   ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
   ! initial perturbation on an isothermal state.
   use grid_data, only: ug,vg,virtempg
   real(r_kind), dimension(nlons,nlats) :: rnh,xnh,rsh,xsh
   real(r_kind) :: lonc,latc,up,pertrad
   integer k
   lonc = pi/9.
   up = 1.
   pertrad = rerth/10.
   latc = 2.*pi/9.
   xnh = sin(latc)*sin(lats) + cos(latc)*cos(lats)*cos(lons-lonc)
   latc = -2.*pi/9.
   xsh = sin(latc)*sin(lats) + cos(latc)*cos(lats)*cos(lons-lonc)
   rnh = rerth*acos(xnh)
   rsh = rerth*acos(xsh)
   virtempg = 300.  ! isothermal state.
   vg = 0.
   do k=1,nlevs
      ! add a barotropic zonal wind perturbation (opp sign in each hemisphere)
      ug(:,:,k) = up*(exp(-(rnh/pertrad)**2)-exp(-(rsh/pertrad)**2))
      call getvrtdivspec(ug(:,:,k),vg(:,:,k),vrtspec(:,k),divspec(:,k),rerth)
      call grdtospec(virtempg(:,:,k),virtempspec(:,k))
   enddo
   topospec = 0.
 end subroutine heldsuarez_ics

end module dyn_init
