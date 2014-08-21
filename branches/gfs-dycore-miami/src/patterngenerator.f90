module patterngenerator

 ! generate random patterns with specified temporal and spatial auto-correlation
 ! in spherical harmonic space.

 use kinds, only: r_kind,r_double,r_single
 use physcons, only:  pi => con_pi, rerth => con_rerth
 use mersenne_twister, only: random_setseed,random_gauss,random_stat
 implicit none
 private

 public :: computevarspec, computevargrid, gaussian_spect, &
  patterngenerator_init, patterngenerator_destroy, getnoise, &
  patterngenerator_advance, getvarspectrum, random_pattern

 type random_pattern
    real(r_kind), public :: lengthscale
    real(r_kind), public :: tau
    real(r_kind), public :: dt
    real(r_kind), public :: phi
    real(r_kind), public :: stdev
    real(r_kind), allocatable, dimension(:), public :: varspectrum, varspectrum1d, lap
    integer, allocatable, dimension(:), public :: degree,order
    integer, public :: seed
    type(random_stat), public :: rstate
 end type random_pattern

 real(r_kind) :: normfact
 integer :: nlons,nlats,ntrunc,ndimspec
  
 contains

 subroutine patterngenerator_init(lscale, delt, tscale, stdev, iseed, rpattern,&
                                  nlon, nlat, jcap, normalization)
   real(r_kind), intent(in) :: lscale,delt,tscale,stdev,normalization
   integer, intent(in) :: nlon,nlat,jcap
   type(random_pattern), intent(out) :: rpattern
   integer, intent(in) :: iseed
   integer m,n
   integer count, count_rate, count_max
   nlons = nlon; nlats = nlat; ntrunc = jcap
   ndimspec = (ntrunc+1)*(ntrunc+2)/2   
   normfact = normalization
   allocate(rpattern%degree(ndimspec),rpattern%order(ndimspec),rpattern%lap(ndimspec))
   rpattern%degree = (/((n,n=m,ntrunc),m=0,ntrunc)/)
   rpattern%order = (/((m,n=m,ntrunc),m=0,ntrunc)/)
   rpattern%lap = -rpattern%degree*(rpattern%degree+1.0)
   rpattern%tau = tscale; rpattern%lengthscale = lscale
   rpattern%dt = delt; rpattern%phi = exp(-delt/tscale)
   rpattern%tau = tscale; rpattern%seed = iseed
   rpattern%stdev = stdev
   allocate(rpattern%varspectrum(ndimspec))
   allocate(rpattern%varspectrum1d(0:ntrunc))
   if (iseed == 0) then
     ! generate a random seed from system clock and ens member number
     call system_clock(count, count_rate, count_max)
     print *,'using seed',count
   else
     count = iseed
     print *,'using seed',count
   endif
   rpattern%seed = count
   call random_setseed(rpattern%seed,rpattern%rstate)
   call gaussian_spect(rpattern)
 end subroutine patterngenerator_init

 subroutine patterngenerator_destroy(rpattern)
   type(random_pattern), intent(inout) :: rpattern
   deallocate(rpattern%varspectrum,rpattern%varspectrum1d)
   deallocate(rpattern%degree,rpattern%order,rpattern%lap)
 end subroutine patterngenerator_destroy

 subroutine computevarspec(rpattern,dataspec,var)
    ! compute globally integrated variance from spectral coefficients
    complex(r_kind), intent(in) :: dataspec(ndimspec)
    real(r_kind), intent(out) ::  var
    type(random_pattern), intent(in) :: rpattern
    integer n
    var = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
           var = var + dataspec(n)*conjg(dataspec(n))
       else
           var = var + 0.5*dataspec(n)*conjg(dataspec(n))
       endif
    enddo
    var = var/normfact
 end subroutine computevarspec

 subroutine computevargrid(datagrid,var,areawts)
    ! compute globally integrated variance from data on gaussian grid
    real(r_kind), intent(in) :: datagrid(nlons,nlats),areawts(nlons,nlats)
    real(r_kind), intent(out) ::  var
    var = sum(areawts*datagrid**2)
 end subroutine computevargrid

 subroutine getvarspectrum(rpattern,dataspec,varspect)
    type(random_pattern), intent(in) :: rpattern
    complex(r_kind), intent(in) :: dataspec(ndimspec)
    real(r_kind), intent(out) :: varspect(0:ntrunc)
    integer n
    varspect = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
          varspect(rpattern%degree(n)) = varspect(rpattern%degree(n)) + &
          dataspec(n)*conjg(dataspec(n))
       else
          varspect(rpattern%degree(n)) = varspect(rpattern%degree(n)) + &
          0.5*dataspec(n)*conjg(dataspec(n))
       endif
    enddo
    varspect = varspect/normfact
 end subroutine getvarspectrum

 subroutine getnoise(rpattern,noise)
   ! generate white noise with unit variance in spectral space
   type(random_pattern), intent(inout) :: rpattern
   complex(r_kind), intent(inout) :: noise(ndimspec)
   real(r_kind) :: noiser(2*ndimspec)
   integer n
   call random_gauss(noiser,rpattern%rstate)
   noiser(1) = 0.; noiser(ndimspec+1) = 0.
   do n=1,ndimspec
      if (rpattern%order(n) .ne. 0.) then
        noise(n) = cmplx(noiser(n),noiser(n+ndimspec))/sqrt(2.*rpattern%degree(n)+1)
      else
        noise(n) = sqrt(2.)*noiser(n)/sqrt(2.*rpattern%degree(n)+1.)
      endif
   enddo
   ! normalize so global mean variance is 1.
   noise = noise*sqrt(normfact/ntrunc)
 end subroutine getnoise

 subroutine patterngenerator_advance(dataspec,rpattern)
    ! advance 1st-order autoregressive process with
    ! specified autocorrelation (phi) and variance spectrum (spectrum)
    complex(r_kind), intent(inout) :: dataspec(ndimspec)
    complex(r_kind) :: noise(ndimspec)
    type(random_pattern), intent(inout) :: rpattern
    call getnoise(rpattern,noise)
    dataspec =  rpattern%phi*dataspec + &
    rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum*noise
 end subroutine patterngenerator_advance

 subroutine gaussian_spect(rpattern)
 ! define variance spectrum (isotropic gaussian covariance)
 ! normalized to unit global variance
  type(random_pattern), intent(inout) :: rpattern
  complex(r_kind) noise(ndimspec)
  real(r_kind) var
  integer n

  ! 1d variance spectrum (as a function of total wavenumber)
  do n=0,ntrunc
     rpattern%varspectrum1d(n) = exp(-rpattern%lengthscale**2*(float(n)*(float(n)+1.))/(4.*rerth**2))
  enddo
  ! scaling factors for spectral coeffs of white noise pattern with unit variance
  rpattern%varspectrum = sqrt(ntrunc*exp(rpattern%lengthscale**2*rpattern%lap/(4.*rerth**2)))
  noise = 0.
  do n=1,ndimspec
     if (rpattern%order(n) .ne. 0.) then
       noise(n) = cmplx(1.,1.)/sqrt(2.*rpattern%degree(n)+1)
     else
       noise(n) = sqrt(2.)/sqrt(2.*rpattern%degree(n)+1.)
     endif
  enddo
  noise(1) = 0 ! no global mean.
  ! normalize so global mean variance is 1.
  noise = noise*sqrt(normfact/ntrunc)
  noise = rpattern%varspectrum*noise
  call computevarspec(rpattern,noise,var)

  ! normalize so patterns will have unit variance
  rpattern%varspectrum = rpattern%varspectrum/sqrt(var)
  rpattern%varspectrum1d = rpattern%varspectrum1d/var

 end subroutine gaussian_spect

end module patterngenerator
