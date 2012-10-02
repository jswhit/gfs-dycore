program testshtns
! ifort -O3 -xHOST -openmp -convert big_endian -traceback testshtns.f90 kinds.o
! shtns.o sigio_module.o ../../lib/libshtns-icc.a ../../intel/lib/libfftw3.a
 use omp_lib, only: omp_get_num_threads, omp_set_num_threads
 use kinds, only: r_kind, r_single, r_double
 use sigio_module, only: sigio_sclose,sigio_swohdc,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_head,sigio_sropen,sigio_srdata,sigio_axdata
 use shtns, only: shtns_init, spectogrd, grdtospec, getgrad, getvrtdivspec, lap, lats, lons, getuv, nlm
 implicit none
 integer, parameter :: ntrunc = 254
 integer, parameter :: ndimspec = (ntrunc+1)*(ntrunc+2)/2
 integer, parameter :: nlons = 768
 integer, parameter :: nlats = 384
 integer, parameter :: ntrac = 3
 integer, parameter :: nlevs = 64
 real(r_double), parameter :: polar_opt=1.d-10
 real(r_kind), parameter :: rerth=6.3712e6
 integer, parameter :: nthreads=1
 character(len=500) :: filename
 type(sigio_data) sigdata
 type(sigio_head) sighead
 integer lu,iret,k,nt,nth
 complex(r_kind), dimension(ndimspec,nlevs) :: vrtspec,divspec,virtempspec
 complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: tracerspec
 complex(r_kind), dimension(ndimspec) ::  lnpsspec,topospec
 real(r_kind), dimension(nlons,nlats,nlevs) :: &
 ug,vg,virtempg,vrtg,divg,dvirtempdx,dvirtempdy
 real(r_kind), dimension(nlons,nlats,nlevs,ntrac) :: tracerg
 real(r_kind), dimension(nlons,nlats) :: lnpsg,topog
 real(8) t1,t2,t0
 integer(8) count, count_rate, count_max

 call getarg(1,filename)
 ! initialize spherical harmonic lib
 call shtns_init(nlons,nlats,ntrunc,nthreads=nthreads,polar_opt=polar_opt)
 ! read spectral initial conditions
 lu = 7
 print *,trim(filename)
 call sigio_srohdc(lu,trim(filename),sighead,sigdata,iret)
 print *,sighead%jcap,sighead%levs
 if (iret .ne. 0) then
   print *,'error reading ',trim(filename),iret
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

 !call OMP_SET_NESTED(.true.)
 do nth=1,24
 call omp_set_num_threads(nth)
 !$omp parallel
 k = omp_get_num_threads()
 !$omp end parallel
 print *,k,'threads'
 call system_clock(count, count_rate, count_max)
 t1 = count*1.d0/count_rate
!$omp parallel do private(k,nt) schedule(dynamic)
 do k=1,nlevs
    call getuv(vrtspec(:,k),divspec(:,k),ug(:,:,k),vg(:,:,k),rerth)
    call spectogrd(vrtspec(:,k),vrtg(:,:,k))
    call spectogrd(divspec(:,k),divg(:,:,k))
    call spectogrd(virtempspec(:,k),virtempg(:,:,k))
    ! gradient of virtual temperature on grid.
    call getgrad(virtempspec(:,k),dvirtempdx(:,:,k),dvirtempdy(:,:,k),rerth)
    ! specific humidity, other tracers on grid.
    do nt=1,ntrac
       call spectogrd(tracerspec(:,k,nt),tracerg(:,:,k,nt))
    enddo
 enddo
!$omp end parallel do 
 call system_clock(count, count_rate, count_max)
 t2 = count*1.d0/count_rate
 if (nth .eq. 1) t0=t2-t1
 print *,'time=',t2-t1,' speedup=',t0/(t2-t1)
 print *,'min/max u,tv:',minval(ug),maxval(ug),minval(virtempg),maxval(virtempg)
 enddo

end program testshtns

subroutine copyspecin(rspecdata,cspecdata,ndimspec)
 use kinds, only: r_kind, r_single, r_double
 implicit none
 integer, intent(in) :: ndimspec
 real(r_single), intent(in) :: rspecdata(2*ndimspec)
 complex(r_kind), intent(out) :: cspecdata(ndimspec)
 real(r_kind) pi
 integer n,nn
 nn = 1
 pi = 4.*atan(1.0)
 ! factor of sqrt(2.*pi) accounts for difference in normalization
 ! between ncep libs and shtns (which uses orthonormalized norm)
 do n=1,ndimspec
    cspecdata(n) = sqrt(2.*pi)*cmplx(rspecdata(nn),rspecdata(nn+1))
    nn = nn + 2
 enddo
end subroutine copyspecin
