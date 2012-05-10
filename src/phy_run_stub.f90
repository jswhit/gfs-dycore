 module phy_run
! do-nothing stub. Provides a dummy getphytend for
! adiabatic model test.

 use params, only: ndimspec, nlevs
 use kinds, only: r_kind
 implicit none
 private

 public :: getphytend

 contains

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt,dlnpsspecdt,dt)
   ! compute physics tendencies for held-suarez test case.
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt
   real(r_kind), intent(in) :: dt
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   return
 end subroutine getphytend

end module phy_run
