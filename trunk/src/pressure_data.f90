module pressure_data
! data for model level pressures and pressure-related
! quantities.
! init_pressdata: allocates arrays.
! calc_pressdata: computes pressure related vars given ln(psg).
! (called by getdyntend in module dyn_run).
! destroy_pressdata: deallocate arrays.
 use params, only: nlons,nlats,nlevs
 use kinds, only: r_kind
 use physcons, only: con_rd,con_cp
 implicit none
 private
 public :: ak,bk,ck,dbk,psg,pk,alfa,rlnp,dpk,&
  prs,init_pressdata, calc_pressdata, destroy_pressdata
! ak,bk: definition of hybrid sigma-pressure coordinates
! dbk(k) = bk(k+1)-bk(k)
! ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
! above 1-d arrays computed in init_dyn, defined from top to bottom
! of model.
! the following are definied by call to calc_pressdata:
! psg: surface pressure (Pa)
! pk: model interface pressures = ak(k) + bk(k)*psg(:,:) (k=1,nlevs+1)
! dpk = pk(:,:,k+1)-pk(:,:,k)
! prs: model layer pressures 
! rlnp = log(pk(:,:,k+1)/pk(:,:,k))
! alfa = 1.-(pk(:,:,k)/dpk(:,:,k))*rlnp(:,:,k) (for k=1 alfa(:,:,1)=log(2))
! all of the above arrays are defined top to bottom, except
! prs which is bottom to top.
 real(r_kind), dimension(:), allocatable :: ak,bk,ck,dbk
 real(r_kind), dimension(:,:), allocatable :: psg
 real(r_kind), dimension(:,:,:), allocatable :: pk,alfa,rlnp,dpk,prs
 
 contains

 subroutine init_pressdata()
    allocate(ak(nlevs+1))
    allocate(bk(nlevs+1))
    allocate(ck(nlevs))
    allocate(dbk(nlevs))
    allocate(psg(nlons,nlats))
    allocate(alfa(nlons,nlats,nlevs))
    allocate(rlnp(nlons,nlats,nlevs))
    allocate(dpk(nlons,nlats,nlevs))
    allocate(pk(nlons,nlats,nlevs+1))
    allocate(prs(nlons,nlats,nlevs))
 end subroutine init_pressdata

 subroutine destroy_pressdata()
    deallocate(ak,bk,ck,dbk)
    deallocate(psg,alfa,rlnp,dpk,pk,prs)
 end subroutine destroy_pressdata

 subroutine calc_pressdata(lnpsg)
    real(r_kind),  intent(in) :: lnpsg(nlons,nlats)
    real(r_kind) rk
! update pressure related variables using latest estimate of lnps
    integer k
    rk = con_rd/con_cp
    ! surface press from log(ps)
    psg = exp(lnpsg)
    ! interface pressure (psg is ps)
    ! ak,bk go top to bottom, so does pk
!$omp parallel
!$omp do private(k)
    do k=1,nlevs+1
      pk(:,:,k)=ak(k) + bk(k)*psg(:,:)
    enddo
!$omp end do
!$omp do private(k)
    do k=1,nlevs
       ! layer pressure thickness
       dpk(:,:,k)=    pk(:,:,k+1) - pk(:,:,k)
       ! sela's layer pressure from hyb2press.f
       ! (goes from bottom to top, unlike pk)
       prs(:,:,nlevs-k+1) = ((pk(:,:,k+1)**rk*pk(:,:,k+1) - pk(:,:,k)**rk*pk(:,:,k))/&
                    ((rk+1.)*dpk(:,:,k))) ** (1./rk)
    enddo
!$omp end do
!$omp end parallel  
    alfa(:,:,1)=log(2.)
    rlnp(:,:,1)=99999.99 !doesn't matter, should never be used.
!$omp parallel do private(k)
    do k=2,nlevs
      rlnp(:,:,k)= log( pk(:,:,k+1)/pk(:,:,k) )
      alfa(:,:,k)= 1.-( pk(:,:,k)/dpk(:,:,k) )*rlnp(:,:,k)
    enddo
!$omp end parallel do 
 end subroutine calc_pressdata

end module pressure_data
