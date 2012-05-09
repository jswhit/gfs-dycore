 module phy_run
! compute physics tendencies for simplifed held-suarez forcing
! getphytend: compute physics tendencies.
! Public subroutines:
! getphytend: compute newtonian damping of temperature, drag terms
! in vort and div eqns.

 use params, only: nlevs,nlons,nlats,ndimspec
 use kinds, only: r_kind
 use shtns, only: grdtospec, lats
 use grid_data, only: vrtg,divg,virtempg
 use pressure_data, only:  prs,psg,pk
 use physcons, only: rd => con_rd, cp => con_cp

 implicit none
 private

 public :: getphytend

 contains

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt,dlnpsspecdt)
   ! compute physics tendencies for held-suarez test case.
   ! http://dx.doi.org/10.1175/1520-0477(1994)075%3C1825:APFTIO%3E2.0.CO;2
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   real(r_kind) p0,kappa,sigbot,tempstrat,delthz,deltmp,&
                kdrag,krada,kradb
   real(r_kind), dimension(:,:,:),allocatable :: blprof,radequiltemp,forcingg
   complex(r_kind), dimension(:,:),allocatable :: forcingspec
   integer k

   allocate(blprof(nlons,nlats,nlevs),radequiltemp(nlons,nlats,nlevs))
   allocate(forcingg(nlons,nlats,nlevs),forcingspec(ndimspec,nlevs))

   kappa = rd/cp
   sigbot = 0.7
   delthz = 10.
   tempstrat = 210.
   kdrag = 1./(1.*86400.)
   krada = 1./(40.*86400.)
   kradb = 1./(4.*86400. )
   p0 = 1.e5
   deltmp = 60.
!$omp parallel do private(k)
   do k=1,nlevs
      blprof(:,:,k) = prs(:,:,k)/psg
      radequiltemp(:,:,k) = (prs(:,:,k)/p0)**(kappa)*&
                     (315.-deltmp*sin(lats)**2-delthz*log(prs(:,:,k)/p0)*(1.+cos(lats)**2))
   enddo
!$omp end parallel do 
   blprof = (blprof-sigbot)/(1.-sigbot)
   where (blprof < 0) 
     blprof = 0
   endwhere
   where (radequiltemp < tempstrat)
      radequiltemp=tempstrat
   endwhere
!$omp parallel do private(k)
   do k=1,nlevs
      forcingg(:,:,k)=(krada+(kradb-krada)*blprof(:,:,k)*cos(lats)**4)*&
                      (radequiltemp(:,:,k)-virtempg(:,:,k))
      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
      dvirtempspecdt(:,k) = dvirtempspecdt(:,k) + forcingspec(:,k)
      forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*vrtg(:,:,k)
      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
      dvrtspecdt(:,k) = dvrtspecdt(:,k) + forcingspec(:,k)
      forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*divg(:,:,k)
      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
      ddivspecdt(:,k) = ddivspecdt(:,k) + forcingspec(:,k)
   enddo
!$omp end parallel do 

   deallocate(blprof,radequiltemp)
   deallocate(forcingg,forcingspec)

   return
 end subroutine getphytend

end module phy_run
