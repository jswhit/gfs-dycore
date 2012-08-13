module finalize_mod
! finalization module (provides finalize driver subroutine)
private
public :: finalize
contains

subroutine finalize
  use dyn_finalize ,only: finalize_dyn
  use phy_finalize ,only: finalize_phy

  implicit none

  call finalize_dyn ()
  call finalize_phy ()
  
  return
end subroutine finalize
end module finalize_mod
