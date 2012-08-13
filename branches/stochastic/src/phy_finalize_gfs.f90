module phy_finalize
! deallocate storage

 use phy_data, only: destroy_phydata
 implicit none
 private
 public :: finalize_phy

 contains

 subroutine finalize_phy()
   call destroy_phydata()
 end subroutine finalize_phy


end module phy_finalize
