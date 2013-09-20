module grid_data
! gaussian grid model state variable data.
! public subroutines:
! init_griddata: allocate arrays.
! destroy_griddata: deallocate arrays.
 use kinds, only: r_kind
 use params, only: nlons,nlats,nlevs,ntrac
 implicit none
 private
 public :: init_griddata, destroy_griddata
 real(r_kind), allocatable, public, dimension(:,:,:) :: ug,vg,vrtg,divg,&
 virtempg,etadot,dlnpdtg
 real(r_kind), allocatable, public, dimension(:,:,:,:) :: tracerg
! (nlons,nlons,nlevs) arrays (bottom to top unless otherwise noted)
! they are transformed to the grid from spectral space in subroutine
! getdyntend in module dyn_run.
! ug: zonal wind 
! vg: meridional wind
! vrtg: vorticity
! divg: divergence
! tracerg: tracers (first one is specific humidity)
! virtempg: virtual temperature
! the following last two are computed in subroutine omega from module dyn_run
! (called by dyntend):
! etadot: vertical motion in hybrid sigma-pressure (top to bottom)
! dlnpdtg: d(lnp)/dt = omega/p (compute in subroutine omega from module dyn_run)
! dlnpsdt: local tendency of ln(ps)
 real(r_kind), allocatable, public, dimension(:,:) :: &
 dlnpsdt,lnpsg,dlnpdst,dphisdx,dphisdy,phis

 contains
 subroutine init_griddata()
    allocate(ug(nlons,nlats,nlevs))
    allocate(vg(nlons,nlats,nlevs))
    allocate(vrtg(nlons,nlats,nlevs))
    allocate(divg(nlons,nlats,nlevs))
    allocate(virtempg(nlons,nlats,nlevs))
    allocate(tracerg(nlons,nlats,nlevs,ntrac))
    allocate(dlnpdtg(nlons,nlats,nlevs))
    allocate(etadot(nlons,nlats,nlevs+1))
    allocate(lnpsg(nlons,nlats))
    allocate(phis(nlons,nlats))
    allocate(dphisdy(nlons,nlats))
    allocate(dphisdx(nlons,nlats))
    allocate(dlnpsdt(nlons,nlats))
 end subroutine init_griddata
 subroutine destroy_griddata()
    deallocate(ug,vg,vrtg,divg,virtempg,dlnpdtg,etadot)
    deallocate(tracerg)
    deallocate(lnpsg,phis,dphisdx,dphisdy,dlnpsdt)
 end subroutine destroy_griddata

end module grid_data
