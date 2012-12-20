module semimp_data
! static data arrays used in semi-implicit time integration scheme.
! Public subroutines:
! init_semimpdata: allocate and populate arrays.
! destroy_semimpdata: deallocate arrays.
 use kinds, only: r_kind, default_real, r_double
 use params, only: dt,ntrunc,nlons,nlats,nlevs
 use pressure_data, only: ak,bk
 use physcons, only: rd => con_rd, cp => con_cp, rerth => con_rerth,&
                     kappa => con_rocp
 implicit none
 private
 public :: init_semimpdata, destroy_semimpdata
 real(r_kind), allocatable, public, dimension(:,:,:) :: d_hyb_m
 real(r_kind), allocatable, public, dimension(:,:) :: amhyb,bmhyb
 real(r_kind), allocatable, public, dimension(:) ::  &
 tref,pkref,dpkref,alfaref,svhyb,tor_hyb
 ! scheme derived from Ascher, Spiteri and Ruuth 1997.
 ! Applied Numerical Mathematics DOI:10.1016/S0168-9274(97)00056-1
 real(r_kind),parameter :: alpha=1.-sqrt(2.)/2.
 ! this choice of delta gives absolute stability for courant numbers < sqrt(2)
 ! it's a combination of ARS 2.5 and 2.6 (tableaus from 2.5, value of delta from 2.6)
 real(r_kind),parameter :: delta=1.-0.5/alpha 
 !real(r_kind),parameter :: delta=-2.*sqrt(2.)/3. ! gives ARS 2.5
 real(r_kind),parameter, public :: a21=alpha
 real(r_kind),parameter, public :: a31=delta
 real(r_kind),parameter, public :: a32=1-delta
 ! Note: it is assumed scheme is diagonally implicit (a22=a33=alpha)
 real(r_kind),parameter, public :: aa21=0.
 real(r_kind),parameter, public :: aa22=alpha
 real(r_kind),parameter, public :: aa31=0.
 real(r_kind),parameter, public :: aa32=1.-alpha
 real(r_kind),parameter, public :: aa33=alpha
 real(r_kind),parameter, public :: b1=0.
 real(r_kind),parameter, public :: b2=1.-alpha
 real(r_kind),parameter, public :: b3=alpha
 real(r_kind),parameter, public :: bb1=0.
 real(r_kind),parameter, public :: bb2=1.-alpha
 real(r_kind),parameter, public :: bb3=alpha
 real(r_kind),parameter, public :: ref_temp = 300.
 real(r_kind),parameter, public :: ref_press = 800.e2

 contains

 subroutine init_semimpdata()
   real(r_double) rnn1
   integer irow,icol,icolbeg,i,j,k,n,icolend,nn,iret
   real(r_double), allocatable, dimension(:,:) :: yecm,tecm,ym,rim
   real(r_double), allocatable, dimension(:) :: vecm
   integer, allocatable, dimension(:) :: ipiv

   ! module vars
   allocate(ipiv(nlevs))
   allocate(amhyb(nlevs,nlevs),bmhyb(nlevs,nlevs))
   allocate(tref(nlevs),pkref(nlevs+1),dpkref(nlevs),alfaref(nlevs))
   allocate(svhyb(nlevs),tor_hyb(nlevs))
   allocate(d_hyb_m(nlevs,nlevs,ntrunc+1))
   ! temp storage.
   allocate(yecm(nlevs,nlevs),tecm(nlevs,nlevs))
   allocate(ym(nlevs,nlevs))
   allocate(vecm(nlevs))
   allocate(rim(nlevs,nlevs)) ! identity matrix

   tref = ref_temp
   do k=1,nlevs+1
      pkref(k) = ak(k) + bk(k)*ref_press
   enddo
   do k=1,nlevs
      dpkref(k) = pkref(k+1)-pkref(k)
   enddo
   alfaref(1) = log(2.)
   do k=2,nlevs
      alfaref(k) = 1.-(pkref(k)/dpkref(k))*log(pkref(k+1)/pkref(k))
   enddo
   yecm=0.
   do irow=1,nlevs
      yecm(irow,irow)=alfaref(irow)*rd
      icolbeg=irow+1
      if(icolbeg<=nlevs)then
       do icol=icolbeg,nlevs
        yecm(irow,icol)=rd*log( pkref(icol+1)/pkref(icol) )
       enddo
      endif
   enddo
   tecm=0.
   do irow=1,nlevs
      tecm(irow,irow)=kappa*tref(irow)*alfaref(irow)
      icolend=irow-1
      do icol=1,icolend
         tecm(irow,icol)=(kappa*tref(irow)*dpkref(icol)/ &
                          dpkref(irow))*log(pkref(irow+1)/pkref(irow))
      enddo
   enddo
   vecm=dpkref/ref_press
   ! amhyb is operator for linearized pres grad term (geopot part)
   ! (operatores on virt temp)
   ! bmhyb is operator for linearized energy conv term
   ! (operates on div)
   do j=1,nlevs
   svhyb(j)=vecm(nlevs+1-j) ! for linearized lnps tend term (div)
   do k=1,nlevs
     amhyb(k,j)=yecm(nlevs+1-k,nlevs+1-j)
     bmhyb(k,j)=tecm(nlevs+1-k,nlevs+1-j)
   enddo
   enddo
   amhyb=amhyb/rerth**2 
   ! times lnps in linearized pgf term (lapacian in div eqn)
   tor_hyb=rd*tref/rerth**2

   rim = 0.
   do k=1,nlevs
      rim(k,k) = 1.
   enddo
! computations that do not depend on wavenumber
   do i=1,nlevs
      do j=1,nlevs
         ym(i,j) = tor_hyb(i)*svhyb(j)
      enddo
      do k=1,nlevs
      do j=1,nlevs
         ym(i,j) = ym(i,j) + amhyb(i,k)*bmhyb(k,j)
      enddo
      enddo
   enddo

! computations that do depend on wavenumber
! enabling openmp for this loop doesn't work with intel MKL
!!$omp parallel do private(nn,n,rnn1,yecm,ipiv,iret,vecm)
   do nn=1,ntrunc+1
      n = nn-1
      rnn1 = n*(n+1)
      yecm = rim + (alpha*dt)**2*rnn1*ym
      ! invert matrix using LAPACK, save in d_hyb_m
      call dgetrf(nlevs,nlevs,yecm,nlevs,ipiv,iret)
      call dgetri(nlevs,yecm,nlevs,ipiv,vecm,nlevs,iret)
      d_hyb_m(:,:,nn) = yecm
   enddo
!!$omp end parallel do 
   deallocate(rim,yecm,tecm,ym,vecm,ipiv)

 end subroutine init_semimpdata

 subroutine destroy_semimpdata()
   deallocate(amhyb,bmhyb,d_hyb_m)
   deallocate(tref,pkref,dpkref,alfaref,svhyb,tor_hyb)
 end subroutine destroy_semimpdata

end module semimp_data
