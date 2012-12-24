module phy_data
! stub
! public subroutines:
! init_phydata: allocate and populate arrays.
! destroy_phydata: deallocate arrays.
 use kinds, only: r_kind,r_double
 use params, only: nlats,nlons
 implicit none
 private
 public :: init_phydata, destroy_phydata, wrtout_flx, wrtout_sfc, pwat
 real(r_kind), allocatable, dimension(:,:) :: pwat

 contains

 subroutine init_phydata()
    allocate(pwat(nlons,nlats))
    pwat = 0.
 end subroutine init_phydata

 subroutine destroy_phydata()
    deallocate(pwat)
 end subroutine destroy_phydata

 subroutine wrtout_sfc(fhour,filename)
   implicit none
   real(r_kind), intent(in) :: fhour
   character(len=120) filename
 end subroutine wrtout_sfc

 subroutine wrtout_flx(fhour,ta,filename)
   implicit none
   real(r_kind),intent(in) :: fhour
   real(r_double),intent(in) :: ta
   character(len=120),intent(in) :: filename
 end subroutine wrtout_flx

end module phy_data
