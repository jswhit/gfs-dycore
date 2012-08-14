program test_patterngenerator
 use kinds, only: r_kind, r_single, r_double
 use shtns, only: shtns_init, spectogrd, grdtospec, getgrad, getvrtdivspec,&
 lap, lats, lons, getuv, nlm, invlap, gauwts, degree, order
 use patterngenerator, only: computevarspec, getnoise, rnorm,&
 patterngenerator_advance, computevargrid, set_random_seed, getvarspectrum,&
 patterngenerator_init,patterngenerator_destroy,random_pattern
 implicit none
 integer, parameter :: nlons=384
 integer, parameter :: nlats=192
 integer, parameter :: ntrunc=126
 real(r_kind), parameter :: rerth=6.3712e6
 real(r_kind), parameter :: dt=1200.
 real(r_kind), parameter :: tau=6.*3600.
 real(r_kind), parameter :: lengthscale=500.e3
 real(r_kind), parameter :: psistdev=1.
 integer, parameter :: nmax=10000
 type(random_pattern) :: rpattern
 real(r_kind) var1,var2,pi
 integer n,iseed,nout
 complex(r_kind), allocatable, dimension(:) :: psispec, noise
 real(r_kind), allocatable, dimension(:,:) :: psigrd,vargrd,&
 psigrd_out,psigrd_save,cov
 real(r_kind) varspect(0:ntrunc),varspect1(0:ntrunc)

 pi = 4.*atan(1.0)

 call shtns_init(nlons,nlats,ntrunc)

 allocate(psispec(nlm),noise(nlm))
 allocate(psigrd(nlons,nlats))
 allocate(psigrd_save(nlons,nlats))
 allocate(psigrd_out(nlons,nlats))
 allocate(vargrd(nlons,nlats))
 allocate(cov(nlons,nlats))

 iseed = 0
 call patterngenerator_init(lengthscale, dt, tau, psistdev, iseed, rpattern)
 print *,'phi = ',rpattern%phi

 ! check spatial variance of red noise random patterns
 call getnoise(psispec)
 psispec = rpattern%varspectrum*psispec
 var2 = 0.
 varspect = 0.
 vargrd = 0.; cov = 0.
 open(7,form='unformatted',access='direct',file='patterns.dat',recl=nlons*nlats*4)
 nout = 1
 do n=1,nmax+1000
    if (n .eq. 999)  call spectogrd(psispec,psigrd)
    call patterngenerator_advance(psispec,rpattern)
    call computevarspec(psispec,var1)
    !print *,n,var1
    if (n .gt. 1000) then
       var2 = var2 + var1/nmax
       psigrd_save = psigrd
       call spectogrd(psispec,psigrd)
       vargrd = vargrd + psigrd**2/nmax
       cov = cov + psigrd*psigrd_save/nmax
       call getvarspectrum(psispec, varspect1)
       varspect = varspect + varspect1/nmax
       ! logit transform
       psigrd_out = (2./(1.+exp(psigrd)))-1.
       nout = nout + 1
       write(7,rec=nout) psigrd_out
    endif
 enddo
 close(7)
 print *,'red noise variance = ',psistdev**2,var2,sum(varspect)
 do n=0,ntrunc
    var1 = rpattern%varspectrum1d(n)*psistdev**2
    print *,n,varspect(n),var1,var1/varspect(n)
 enddo
 print *,'min/max/mean vargrd =',&
 minval(vargrd),maxval(vargrd),sum(rpattern%areawts*vargrd)
 cov = cov/vargrd
 print *,'min/max/mean phi =',&
 minval(cov),maxval(cov),sum(rpattern%areawts*cov)
 
 deallocate(psispec,noise,psigrd,vargrd,psigrd_save,psigrd_out,cov)

 call patterngenerator_destroy(rpattern)

end program test_patterngenerator
