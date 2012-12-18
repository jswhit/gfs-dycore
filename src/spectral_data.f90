module spectral_data
! spectral model state variable data
! Public subroutines:
! init_specdata: allocate arrays.
! destroy_specdata: deallocate arrays.
 use kinds, only: r_kind
 use params, only: ndimspec, nlevs, ntrac
 implicit none
 private
 public :: init_specdata,destroy_specdata
 complex(r_kind), allocatable, public, dimension(:,:) :: &
 vrtspec,divspec,virtempspec
 complex(r_kind), allocatable, public, dimension(:,:,:) :: tracerspec
! ndimspec by nlevs complex arrays (ndimspec = (ntrunc+1)*(ntrunc)/2):
! these arrays are update each timestep in subroutine advance in run.f90.
! vrtspec: vorticity
! divspec: divergence
! virtempspec: virtural temp
! tracerspec: tracers (ntrac=1 is  specific humidity, ntrac=ntoz ozone, ntrac=ntclw cloud
! condensate)
 complex(r_kind), allocatable, public, dimension(:) :: lnpsspec,topospec
 real(r_kind), allocatable, public, dimension(:) :: disspec,diff_prof,dmp_prof
! ndimspec 1-d arrays:
! lnpsspec: log(ps) in Pa
! topospec: orography in meters (static, created in init_dyn)
! disspec: hyperdiffusion operator (static, created in init_dyn)
! diff_prof: vertical profile of hyperdiffusion coefficient.
! dmp_prof: vertical profile of linear momentum drag (active
! at top of model).

 contains
 subroutine init_specdata()
    allocate(vrtspec(ndimspec,nlevs))
    allocate(divspec(ndimspec,nlevs))
    allocate(tracerspec(ndimspec,nlevs,ntrac))
    allocate(virtempspec(ndimspec,nlevs))
    allocate(topospec(ndimspec))
    allocate(lnpsspec(ndimspec))
    allocate(disspec(ndimspec),dmp_prof(nlevs),diff_prof(nlevs))
 end subroutine init_specdata
 subroutine destroy_specdata()
    deallocate(vrtspec,divspec,virtempspec)
    deallocate(tracerspec)
    deallocate(topospec,lnpsspec,disspec,dmp_prof,diff_prof)
 end subroutine destroy_specdata

end module spectral_data
