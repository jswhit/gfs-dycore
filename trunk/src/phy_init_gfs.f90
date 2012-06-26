module phy_init
! stub that provides do nother init_phy routine.

 use phy_data, only: init_phydata
 use funcphys, only: gfuncphys
 implicit none
 private
 public :: init_phy

 contains

 subroutine init_phy()
   call gfuncphys()
   call init_phydata()
 end subroutine init_phy

end module phy_init
