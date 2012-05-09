module spectral_data
! spectral model state variable data
! Public subroutines:
! init_specdata: allocate arrays.
! destroy_specdata: deallocate arrays.
 use kinds, only: r_kind
 use params, only: ndimspec, nlevs
 implicit none
 private
 public :: init_specdata,destroy_specdata
 complex(r_kind), allocatable, public, dimension(:,:) :: &
 vrtspec,divspec,virtempspec,spfhumspec
! ndimspec by nlevs complex arrays (ndimspec = (ntrunc+1)*(ntrunc)/2):
! vrtspec: vorticity
! divspec: divergence
! virtempspec: virtural temp
! spfhumspec: specific humidity
 complex(r_kind), allocatable, public, dimension(:) :: lnpsspec,topospec,&
 disspec
! ndimspec 1-d arrays:
! lnpsspec: log(ps) in Pa
! topospec: orography in meters.
! disspec: hyperdiffusion operator

 contains
 subroutine init_specdata()
    allocate(vrtspec(ndimspec,nlevs))
    allocate(divspec(ndimspec,nlevs))
    allocate(spfhumspec(ndimspec,nlevs))
    allocate(virtempspec(ndimspec,nlevs))
    allocate(topospec(ndimspec))
    allocate(lnpsspec(ndimspec))
    allocate(disspec(ndimspec))
 end subroutine init_specdata
 subroutine destroy_specdata()
    deallocate(vrtspec,divspec,virtempspec,spfhumspec)
    deallocate(topospec,lnpsspec,disspec)
 end subroutine destroy_specdata

end module spectral_data
