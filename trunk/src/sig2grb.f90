program sig2grb

! utility to convert GFS spectral file to grib1 file.

! grib1 file includes ps,zs,p,dp,u,v,tv,q,ozone,clwmr,omega
! on model levels.

! usage:  sig2grb <spectral file> <grib file>

! Jeffrey.S.Whitaker <jeffrey.s.whitaker@noaa.gov> January, 2013

 use kinds, only: r_kind, r_single, r_double
 use sigio_module, only: sigio_sclose,sigio_srhead,&
  sigio_srohdc,sigio_data,sigio_head,sigio_sropen,sigio_axdata
 use shtns, only: shtns_init, spectogrd, grdtospec, getgrad, getvrtdivspec, getuv
 use physcons, only: rerth => con_rerth
 implicit none

 character(len=500) sigfile, gribfile
 type(sigio_data) sigdata
 type(sigio_head) sighead
 integer lu,iret,k,nt,nlons,nlats,nlevs,ntrunc,ntrac,ndimspec
 complex(r_kind), allocatable, dimension(:,:) :: vrtspec,divspec,virtempspec
 complex(r_kind), allocatable, dimension(:,:,:):: tracerspec
 complex(r_kind), allocatable, dimension(:) ::  lnpsspec,topospec
 real(r_kind), allocatable, dimension(:,:,:) :: &
 ug,vg,virtempg,divg
 real(r_kind), allocatable, dimension(:,:,:,:) :: tracerg
 real(r_kind), allocatable, dimension(:,:) :: psg,topog,dlnpsdx,dlnpsdy
 real(8) t1,t2
 integer(8) count, count_rate, count_max

 lu = 7

 call getarg(1,sigfile)
 call getarg(2,gribfile)
 ! read header from spectral file.
 call sigio_sropen(lu,trim(sigfile),iret)
 if (iret .ne. 0) then
    print *,'error opening ',trim(sigfile),iret
    stop
 endif
 call sigio_srhead(lu,sighead,iret)
 if (iret .ne. 0) then
    print *,'error reading header from ',trim(sigfile),iret
    stop
 else
    nlons = sighead%lonb
    nlats = sighead%latb
    nlevs = sighead%levs
    ntrunc = sighead%jcap
    ndimspec = (ntrunc+1)*(ntrunc+2)/2
    ntrac = sighead%ntrac
    print *,'nlons,nlats,nlevs,ntrac,ntrunc=',nlons,nlats,nlevs,ntrac,ntrunc
 endif 
 ! initialize spherical harmonic lib
 call shtns_init(nlons,nlats,ntrunc)
 allocate(vrtspec(ndimspec,nlevs))
 allocate(divspec(ndimspec,nlevs))
 allocate(virtempspec(ndimspec,nlevs))
 allocate(tracerspec(ndimspec,nlevs,ntrac))
 allocate(topospec(ndimspec),lnpsspec(ndimspec))
 allocate(ug(nlons,nlats,nlevs))
 allocate(vg(nlons,nlats,nlevs))
 allocate(divg(nlons,nlats,nlevs))
 allocate(virtempg(nlons,nlats,nlevs))
 allocate(tracerg(nlons,nlats,nlevs,ntrac))
 allocate(psg(nlons,nlats),topog(nlons,nlats),dlnpsdx(nlons,nlats),dlnpsdy(nlons,nlats))
 ! read spectral initial conditions
 lu = 7
 print *,'reading... ',trim(sigfile)
 call sigio_srohdc(lu,trim(sigfile),sighead,sigdata,iret)
 if (iret .ne. 0) then
   print *,'error reading ',trim(sigfile),iret
   stop
 endif
 ! convert spectral arrays to double precision complex,
 ! re-normalize coefficients.
 call copyspecin(sigdata%ps, lnpsspec,ndimspec)
 call copyspecin(sigdata%hs, topospec,ndimspec)
 do k=1,nlevs
    call copyspecin(sigdata%z(:,k),vrtspec(:,k),ndimspec)
    call copyspecin(sigdata%d(:,k),divspec(:,k),ndimspec)
    call copyspecin(sigdata%t(:,k),virtempspec(:,k),ndimspec)
    do nt=1,ntrac
       call copyspecin(sigdata%q(:,k,nt),tracerspec(:,k,nt),ndimspec)
    enddo
 enddo
 call sigio_axdata(sigdata,iret)
 call sigio_sclose(lu,iret)

 call system_clock(count, count_rate, count_max)
 t1 = count*1.d0/count_rate
!$omp parallel do private(k,nt)
 do k=1,nlevs
    call getuv(vrtspec(:,k),divspec(:,k),ug(:,:,k),vg(:,:,k),rerth)
    call spectogrd(divspec(:,k),divg(:,:,k))
    call spectogrd(virtempspec(:,k),virtempg(:,:,k))
    ! specific humidity, other tracers on grid.
    do nt=1,ntrac
       call spectogrd(tracerspec(:,k,nt),tracerg(:,:,k,nt))
    enddo
 enddo
!$omp end parallel do
 call system_clock(count, count_rate, count_max)
 t2 = count*1.d0/count_rate
 print *,'time to do transforms = ',t2-t1
 print *,'min/max tv',minval(virtempg),maxval(virtempg)
 call spectogrd(lnpsspec,psg)
 call spectogrd(topospec,topog)
 psg = 10.*exp(psg) ! convert to Pa.
 print *,'min/max ps',minval(psg),maxval(psg)
 call getgrad(lnpsspec, dlnpsdx, dlnpsdy, rerth)
 print *,'writing.... ',trim(gribfile)
 call wrtout_gfsgrb(sighead,nlons,nlats,nlevs,ntrac,ntrunc,&
                    psg,dlnpsdx,dlnpsdy,topog,ug,vg,divg,virtempg,tracerg,gribfile)

 deallocate(vrtspec,divspec,virtempspec,tracerspec,topospec,lnpsspec)
 deallocate(ug,vg,divg,virtempg,tracerg)
 deallocate(psg,topog,dlnpsdx,dlnpsdy)

 end program sig2grb
 
 subroutine copyspecin(rspecdata,cspecdata,ndimspec)
    use kinds, only: r_single, r_kind
    use physcons, only: pi => con_pi
    implicit none
    integer, intent(in) :: ndimspec
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

 subroutine wrtout_gfsgrb(sighead,nlons,nlats,nlevs,ntrac,ntrunc,&
                          psg,dlnpsdx,dlnpsdy,topog,ug,vg,divg,virtempg,tracerg,filename)
    use kinds, only: r_single, r_kind
    use sigio_module, only: sigio_head
    use gfsio_module, only: gfsio_gfile, gfsio_open, gfsio_writerecvw34, &
                            gfsio_init, gfsio_close
    use physcons, only: con_rd,con_cp,rk => con_rocp, fv => con_fvirt
    implicit none
    type(sigio_head), intent(in) :: sighead
    integer, intent(in) :: nlons,nlats,nlevs,ntrunc,ntrac
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: &
     ug,vg,divg,virtempg
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs,ntrac) :: tracerg
    real(r_kind), intent(in), dimension(nlons,nlats) :: &
     psg,topog,dlnpsdx,dlnpsdy
    ! write out gfsio grib data
    type(gfsio_gfile)    :: gfile
    character(len=500), intent(in) :: filename
    integer k,iret,reclev(2+nlevs*8)
    character(len=8) recname(2+nlevs*8)
    character(len=16) reclevtyp(2+nlevs*9)
    real(4) tmpg(nlons,nlats)
    real(r_kind) prs(nlons,nlats,nlevs),dpk(nlons,nlats,nlevs),pk(nlons,nlats,nlevs+1),&
    ak(nlevs),bk(nlevs),dlnpsdt(nlons,nlats),etadot(nlons,nlats,nlevs+1),&
    dlnpdtg(nlons,nlats,nlevs),dbk(nlevs),ck(nlevs)
    real(8) t1,t2
    integer(8) count, count_rate, count_max

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
    do k=1,nlevs+1
       pk(:,:,k)=ak(k) + bk(k)*psg(:,:)
    enddo
    do k=1,nlevs
       ! layer pressure thickness
       dpk(:,:,k)=    pk(:,:,k+1) - pk(:,:,k)
       ! sela's layer pressure from hyb2press.f
       ! (goes from bottom to top, unlike pk)
       prs(:,:,nlevs-k+1) = ((pk(:,:,k+1)**rk*pk(:,:,k+1) - pk(:,:,k)**rk*pk(:,:,k))/&
                    ((rk+1.)*dpk(:,:,k))) ** (1./rk)
    enddo
    ! compute vertical velocity.
    call system_clock(count, count_rate, count_max)
    t1 = count*1.d0/count_rate
    call getomega(nlons,nlats,nlevs,ug,vg,divg,bk,ck,dbk,pk,dpk,psg,dlnpsdx,dlnpsdy,dlnpsdt,dlnpdtg,etadot)
    call system_clock(count, count_rate, count_max)
    t2 = count*1.d0/count_rate
    print *,'time to compute vertical velocity = ',t2-t1

    recname(1)   = 'hgt'
    reclevtyp(1) = 'sfc'
    reclev(1)    = 1
    recname(2)   = 'pres'
    reclevtyp(2) = 'sfc'
    reclev(2)    = 1
    do k=1,nlevs
      recname(k+2)            = 'pres'
      reclevtyp(k+2)          = 'layer'
      reclev(k+2)             = k
      recname(k+2+nlevs)      = 'dpres'
      reclevtyp(k+2+nlevs)    = 'layer'
      reclev(k+2+nlevs)       = k
      recname(k+2+nlevs*2)    = 'tmp'
      reclevtyp(k+2+nlevs*2)  = 'layer'
      reclev(k+2+nlevs*2)     = k
      recname(k+2+nlevs*3)    = 'ugrd'
      reclevtyp(k+2+nlevs*3)  = 'layer'
      reclev(k+2+nlevs*3)     = k
      recname(k+2+nlevs*4)    = 'vgrd'
      reclevtyp(k+2+nlevs*4)  = 'layer'
      reclev(k+2+nlevs*4)     = k
      recname(k+2+nlevs*5)    = 'spfh'
      reclevtyp(k+2+nlevs*5)  = 'layer'
      reclev(k+2+nlevs*5)     = k
      recname(k+2+nlevs*6)    = 'o3mr'
      reclevtyp(k+2+nlevs*6)  = 'layer'
      reclev(k+2+nlevs*6)     = k
      recname(k+2+nlevs*7)    = 'clwmr'
      reclevtyp(k+2+nlevs*7)  = 'layer'
      reclev(k+2+nlevs*7)     = k
      recname(k+2+nlevs*8)    = 'vvel'
      reclevtyp(k+2+nlevs*8)  = 'layer'
      reclev(k+2+nlevs*8)     = k
    enddo

    !print *,' calling gfsio_open idate=',idate_start,' fhour=',fhour4
 
    call gfsio_init(iret)
    call gfsio_open(gfile,trim(filename),'write',iret,&
         version=sighead%ivs,fhour=sighead%fhour,idate=sighead%idate,nrec=2+nlevs*9,&
         latb=nlats,lonb=nlons,levs=nlevs,jcap=ntrunc,itrun=sighead%itrun,&
         iorder=sighead%iorder,irealf=sighead%irealf,igen=sighead%igen,latf=nlats,lonf=nlons,&
         latr=nlats,lonr=nlons,ntrac=ntrac,icen2=sighead%icen2,iens=sighead%iens,&
         idpp=sighead%idpp,idsl=sighead%idsl,idvc=sighead%idvc,idvm=sighead%idvm,&
         idvt=sighead%idvt,idrun=sighead%idrun,&
         idusr=sighead%idusr,pdryini=sighead%pdryini,ncldt=sighead%ncldt,nvcoord=sighead%nvcoord,&
         vcoord=sighead%vcoord,recname=recname,reclevtyp=reclevtyp,&
         reclev=reclev,Cpi=sighead%cpi,Ri=sighead%ri)

    call twodtooned(topog,tmpg,nlons,nlats)
    call gfsio_writerecvw34(gfile,'hgt','sfc',1,tmpg,iret)
    call twodtooned(psg,tmpg,nlons,nlats)
    call gfsio_writerecvw34(gfile,'pres','sfc',1,tmpg,iret)
    do k=1,nlevs
       call twodtooned(prs(:,:,k),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'pres','layer',k,&
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(dpk(:,:,nlevs-k+1),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'dpres','layer',k, &
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(virtempg(:,:,k),tmpg,nlons,nlats)
       tmpg = tmpg/(1.+fv*tracerg(:,:,k,1))
       call gfsio_writerecvw34(gfile,'tmp','layer',k,&
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(ug(:,:,k),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'ugrd','layer',k, &
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(vg(:,:,k),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'vgrd','layer',k,&
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(tracerg(:,:,k,1),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'spfh','layer',k,&
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(tracerg(:,:,k,2),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'o3mr','layer',k,&
                            tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(tracerg(:,:,k,3),tmpg,nlons,nlats)
       call gfsio_writerecvw34(gfile,'clwmr','layer',k,&
                               tmpg, iret)
    enddo
    do k=1,nlevs
       call twodtooned(dlnpdtg(:,:,k),tmpg,nlons,nlats)
       !tmpg = tmpg*prs(:,:,k)
       tmpg = tmpg*0.5*(pk(:,:,nlevs-k+2)+pk(:,:,nlevs-k+1))
       call gfsio_writerecvw34(gfile,'vvel','layer',k,&
                               tmpg, iret, precision=6)
    enddo

    call gfsio_close(gfile,iret)

 end subroutine wrtout_gfsgrb

 subroutine twodtooned(data2,data1,nlons,nlats)
   use kinds, only: r_kind
   implicit none
   integer, intent(in) :: nlons,nlats
   real(r_kind), intent(in) :: data2(nlons,nlats)
   real(4), intent(out) :: data1(nlons*nlats)
   integer i,j,n
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      data1(n) = data2(i,j)
   enddo
 end subroutine twodtooned

 subroutine getomega(nlons,nlats,nlevs,ug,vg,divg,bk,ck,dbk,pk,dpk,psg,dlnpsdx,dlnpsdy,dlnpsdt,dlnpdtg,etadot)
    use kinds, only: r_kind
    implicit none
    ! compute omega, etadot, tendency of lnps 
    ! all input and output arrays oriented bottom to top (k=1 is near ground)
    integer, intent(in) :: nlons,nlats,nlevs
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: ug,vg,divg,dpk      
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs+1) :: pk
    real(r_kind), intent(in), dimension(nlons,nlats) :: psg,dlnpsdx,dlnpsdy
    real(r_kind), intent(in), dimension(nlevs+1) :: bk
    real(r_kind), intent(in), dimension(nlevs) :: ck,dbk
    ! omega (pressure vertical velocity divided by pressure) on model layers.
    real(r_kind), intent(out), dimension(nlons,nlats,nlevs) :: dlnpdtg
    ! etadot (vertical velocity in hybrid coords) on layer interfaces.
    real(r_kind), intent(out), dimension(nlons,nlats,nlevs+1) :: etadot
    real(r_kind), intent(inout), dimension(nlons,nlats) :: dlnpsdt
! work space:
    real(r_kind), dimension(nlons,nlats,nlevs) :: &
    workb,workc,cg,cb,db,alfa,rlnp
! local scalars 
    integer k

    alfa(:,:,1)=log(2.)
    rlnp(:,:,1)=99999.99 !doesn't matter, should never be used.
    do k=2,nlevs
      rlnp(:,:,k)= log( pk(:,:,k+1)/pk(:,:,k) )
      alfa(:,:,k)= 1.-( pk(:,:,k)/dpk(:,:,k) )*rlnp(:,:,k)
    enddo

!$omp parallel do private(k)
    do k=1,nlevs
     cg(:,:,k)=ug(:,:,nlevs+1-k)*dlnpsdx(:,:)+vg(:,:,nlevs+1-k)*dlnpsdy(:,:)
    enddo
!$omp end parallel do 
 
    db(:,:,1)=divg(:,:,nlevs)*dpk(:,:,1)
    cb(:,:,1)=cg(:,:,1)*dbk(1)
 
    do k=1,nlevs-1
     db(:,:,k+1)=db(:,:,k)+divg(:,:,nlevs-k)*dpk(:,:,k+1)
     cb(:,:,k+1)=cb(:,:,k)+cg(:,:,k+1)*dbk(k+1)
    enddo

    dlnpsdt(:,:) = -db(:,:,nlevs)/psg(:,:)-cb(:,:,nlevs)
    etadot(:,:,1) = 0.; etadot(:,:,nlevs+1) = 0.
!$omp parallel do private(k)
    do k=1,nlevs-1
       etadot(:,:,k+1)=-psg(:,:)*(bk(k+1)*dlnpsdt(:,:)+cb(:,:,k)) - db(:,:,k) 
    enddo
!$omp end parallel do 
 
    workb(:,:,1)=alfa(:,:,1)* &
                ( divg(:,:,nlevs)*dpk(:,:,1)+psg(:,:)*cb(:,:,1)*dbk(1) )
 
!$omp parallel do private(k)
    do k=2,nlevs
      workb(:,:,k)=rlnp(:,:,k)*( db(:,:,k-1)+psg(:,:)*cb(:,:,k-1) )+&
      alfa(:,:,k)*( divg(:,:,nlevs+1-k)*dpk(:,:,k)+psg(:,:)*cg(:,:,k)*dbk(k) )
    enddo
!$omp end parallel do 
 
    workc(:,:,1)=psg(:,:)*cg(:,:,1)*dbk(1)

!$omp parallel
!$omp do private(k)
    do k=2,nlevs
      workc(:,:,k)=psg(:,:)*cg(:,:,k)*( dbk(k)+ck(k)*rlnp(:,:,k)/dpk(:,:,k) )
    enddo
!$omp end do
!$omp do private(k)
    do k=1,nlevs
      dlnpdtg(:,:,nlevs+1-k)=(workc(:,:,k)-workb(:,:,k))/dpk(:,:,k)
    enddo
!$omp end do
!$omp end parallel 
 
    return
 end subroutine getomega 
