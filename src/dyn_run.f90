 module dyn_run
! calculation of model tendencies.
! Public subroutines:
! getdyntend:  dynamics tendencies (not including linear damping and diffusion)
! (called in run.f90)
! semimpadj:  adjustment to tendencies from semi-implicit step.
! (only called in run.f90 if explicit=.false.)
! used by getdyntend:
! getomega: calculate dlnpdt (omega/p)
! getpresgrad: calculate pressure gradient terms in momentum eqns.
! getvadv: calculate vertical advection terms.

 use params, only: nlevs,ntrunc,nlons,nlats,ndimspec,dry,dt,ntrac
 use kinds, only: r_kind,r_double
 use shtns, only: degree,order,&
 lap,invlap,lats,grdtospec,spectogrd,getuv,getvrtdivspec,getgrad
 use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
 tracerspec, topospec
 use grid_data, only: ug,vg,vrtg,virtempg,divg,tracerg,dlnpdtg,etadot,lnpsg,&
 phis,dphisdx,dphisdy,dlnpsdt
 use pressure_data, only:  ak, bk, ck, dbk, dpk, rlnp, pk, alfa, dpk, psg,&
 calc_pressdata
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
               eps => con_eps, omega => con_omega, cvap => con_cvap, &
               grav => con_g, fv => con_fvirt, kappa => con_rocp
 use semimp_data, only: amhyb,bmhyb,svhyb,d_hyb_m,tor_hyb,tref,pkref,alfaref

 implicit none
 private

 public :: getdyntend, getomega, getpresgrad, getvadv, semimpadj

 contains

 subroutine getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,&
                       dtracerspecdt,dlnpsspecdt,just_do_inverse_transform)
   ! compute dynamics tendencies (notincluding hyper-diffusion.and linear drag,
   ! which are treated implicitly)
   ! based on hybrid sigma-pressure dynamical core described in
   ! http://www.emc.ncep.noaa.gov/officenotes/newernotes/on462.pdf
   ! The only difference is that here we use a forward in time
   ! runge-kutta scheme instead of assellin-filtered leap-frog.
   logical, optional, intent(in) :: just_do_inverse_transform
   logical :: early_return = .false. ! if true. spectral-> grid only
   complex(r_kind), intent(out), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
   dtracerspecdt
   complex(r_kind), intent(out), dimension(ndimspec) :: dlnpsspecdt
! local variables.
   complex(r_kind), dimension(:,:),allocatable :: workspec
   real(r_kind), dimension(:,:), allocatable :: dlnpsdx,dlnpsdy
   real(r_kind), dimension(:,:,:), allocatable :: &
   prsgx,prsgy,vadvu,vadvv,vadvt,vadvq,dvirtempdx,dvirtempdy
   integer k,nt
   real(8) t1,t2,t0
   integer(8) count, count_rate, count_max

   ! should tendencies be computed, or just spectral -> grid?
   if (present(just_do_inverse_transform)) then
     early_return = .true.
   else
     early_return = .false.
   endif

   ! use alloctable arrays instead of automatic arrays 
   ! for these to avoid segfaults in openmp.
   allocate(dlnpsdx(nlons,nlats))
   allocate(dlnpsdy(nlons,nlats))
   allocate(prsgx(nlons,nlats,nlevs))
   allocate(prsgy(nlons,nlats,nlevs))
   allocate(vadvu(nlons,nlats,nlevs))
   allocate(vadvv(nlons,nlats,nlevs))
   allocate(vadvt(nlons,nlats,nlevs))
   allocate(vadvq(nlons,nlats,nlevs))
   allocate(dvirtempdx(nlons,nlats,nlevs))
   allocate(dvirtempdy(nlons,nlats,nlevs))
   allocate(workspec(ndimspec,nlevs))  

   ! compute u,v,virt temp, vorticity, divergence, ln(ps)
   ! and specific humidity on grid from spectral coefficients.
   call system_clock(count, count_rate, count_max)
   t1 = count*1.d0/count_rate
   t0=t1
!$omp parallel do private(k)
   do k=1,nlevs
      call getuv(vrtspec(:,k),divspec(:,k),ug(:,:,k),vg(:,:,k),rerth)
      call spectogrd(vrtspec(:,k),vrtg(:,:,k))
      call spectogrd(divspec(:,k),divg(:,:,k))
      call spectogrd(virtempspec(:,k),virtempg(:,:,k))
      ! gradient of virtual temperature on grid.
      call getgrad(virtempspec(:,k),dvirtempdx(:,:,k),dvirtempdy(:,:,k),rerth)
      !print *,k,maxval(abs(ug(:,:,k))),maxval(abs(virtempg(:,:,k))),maxval(abs(dvirtempdx(:,:,k)))
      ! specific humidity, other tracers on grid.
      do nt=1,ntrac
         call spectogrd(tracerspec(:,k,nt),tracerg(:,:,k,nt))
      enddo
   enddo
!$omp end parallel do 
   call system_clock(count, count_rate, count_max)
   t2 = count*1.d0/count_rate
   !print *,'1 time=',t2-t1
   ! lnps on grid.
   call spectogrd(lnpsspec, lnpsg)
   ! gradient of surface pressure.
   ! gradient of surface geopotential precomputed in dyn_init,
   ! saved in module grid_data (dphisdx,dphisdy).
   call getgrad(lnpsspec, dlnpsdx, dlnpsdy, rerth)
   ! compute pressure on interfaces and model layers using lnpsg
   ! i.e. calculate alpha,delp,rlnp,dpk,psg,pk from ak,bk,lnpsg
   ! results stored in module pressure_data
   call calc_pressdata(lnpsg) 
   !print *,'min/max ps',minval(psg/100.),maxval(psg/100.)
   ! get pressure vertical velocity divided by pressure (dlnpdtg),
   ! tendency of lnps, etadot.
   ! etadot goes top to bottom, dlnpdtg goes bottom to top 
   call getomega(ug,vg,divg,dlnpsdx,dlnpsdy,dlnpsdt,dlnpdtg,etadot,&
                 vadvu,vadvv,vadvt,vadvq,prsgx) ! work storage
   ! return before computing tendencies (after updating
   ! grid data).
   if (early_return) then
      deallocate(vadvq,workspec,dvirtempdx,dvirtempdy)
      deallocate(prsgx,prsgy,vadvu,vadvv,vadvt)
      deallocate(dlnpsdx,dlnpsdy)
      return
   endif
   call grdtospec(dlnpsdt,dlnpsspecdt)
   ! get pressure gradient terms (prsgx,prsgy)
   call getpresgrad(virtempg,dvirtempdx,dvirtempdy,dphisdx,dphisdy,dlnpsdx,dlnpsdy,&
                    prsgx,prsgy,&
                    vadvu,vadvv,vadvt,vadvq) ! work storage
   ! get vertical advection terms  (input etadot is top to bottom)
   call getvadv(ug,etadot,vadvu)
   call getvadv(vg,etadot,vadvv)
   call getvadv(virtempg,etadot,vadvt)
   ! add pressure gradient force to vertical advection terms
   ! compute energy conversion term.
   dlnpsdx = 2.*omega*sin(lats) ! temp storage of planetary vorticity
   if (ntrac > 0) then
      ! energy conversion term, store in vadvq
      !$omp parallel do private(k)
      do k=1,nlevs
         ! section 1.5 in ON 461 (eqn 1.0.3).
         !vadvq(:,:,k) = &
         !kappa*dlnpdtg(:,:,k)*virtempg(:,:,k)*&
         !(1.+fv*tracerg(:,:,k,1))/(1.+(cvap/cp-1.)*tracerg(:,:,k,1))
         ! GFS appears to missing term in numerator (1. + fv*q)?
         ! GFS has (in gfidi_hyb.f).  This form is consistent with
         ! ECMWF IFS documentation.
         vadvq(:,:,k) = &
         kappa*dlnpdtg(:,:,k)*virtempg(:,:,k)/(1.+(cvap/cp-1.)*tracerg(:,:,k,1))
      enddo
      !$omp end parallel do 
   else
      !$omp parallel do private(k)
      do k=1,nlevs
         vadvq(:,:,k) = kappa*dlnpdtg(:,:,k)*virtempg(:,:,k)
      enddo
      !$omp end parallel do 
   endif
   call system_clock(count, count_rate, count_max)
   t1 = count*1.d0/count_rate
   !print *,'2 time=',t1-t2
   ! compute tendencies of virt temp, ort and div in spectral space
!$omp parallel do private(k)
   do k=1,nlevs
      ! add pressure gradient force to vertical advection terms
      ! (so prsgy and prsgx can be re-used)
      vadvv(:,:,k) = vadvv(:,:,k) - prsgy(:,:,k)
      vadvu(:,:,k) = vadvu(:,:,k) - prsgx(:,:,k)
      ! virtual temp tendency
      prsgx(:,:,k) = -ug(:,:,k)*dvirtempdx(:,:,k) - vg(:,:,k)*dvirtempdy(:,:,k) - &
                      vadvt(:,:,k) + vadvq(:,:,k)
      call grdtospec(prsgx(:,:,k), dvirtempspecdt(:,k))
      ! flux terms for vort, div eqns
      prsgx(:,:,k) = ug(:,:,k)*(vrtg(:,:,k) + dlnpsdx(:,:)) + vadvv(:,:,k)
      prsgy(:,:,k) = vg(:,:,k)*(vrtg(:,:,k) + dlnpsdx(:,:)) - vadvu(:,:,k)
      call getvrtdivspec(prsgx(:,:,k),prsgy(:,:,k),ddivspecdt(:,k),dvrtspecdt(:,k),rerth)
      ! flip sign of vort tend.
      dvrtspecdt(:,k) = -dvrtspecdt(:,k) 
      ! add laplacian(KE) term to div tendency
      prsgx(:,:,k) = 0.5*(ug(:,:,k)**2+vg(:,:,k)**2)
      call grdtospec(prsgx(:,:,k),workspec(:,k))
      ddivspecdt(:,k) = ddivspecdt(:,k) - &
      (lap(:)/rerth**2)*workspec(:,k)
   enddo
!$omp end parallel do 
   call system_clock(count, count_rate, count_max)
   t2 = count*1.d0/count_rate
   !print *,'3 time=',t2-t1
   ! compute tendency of tracers (including specific humidity) in spectral space.
   do nt=1,ntrac
   ! use positive-definite vertical advection.
   call getvadv_tracers(tracerg(:,:,:,nt),etadot,vadvq)
!$omp parallel do private(k)
   do k=1,nlevs
      ! gradient of specific humidity on grid.
      call getgrad(tracerspec(:,k,nt),dvirtempdx(:,:,k),dvirtempdy(:,:,k),rerth)
      ! specific humidity tendency
      prsgx(:,:,k) = &
      -ug(:,:,k)*dvirtempdx(:,:,k)-vg(:,:,k)*dvirtempdy(:,:,k)-vadvq(:,:,k)
      call grdtospec(prsgx(:,:,k), dtracerspecdt(:,k,nt))
   enddo
!$omp end parallel do 
   enddo
   call system_clock(count, count_rate, count_max)
   t1 = count*1.d0/count_rate

   deallocate(vadvq,workspec,dvirtempdx,dvirtempdy)
   deallocate(prsgx,prsgy,vadvu,vadvv,vadvt)
   deallocate(dlnpsdx,dlnpsdy)

   return
 end subroutine getdyntend

 subroutine semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                      divspec_prev,virtempspec_prev,lnpsspec_prev,kt,dtx)
! semi-implicit adjustment of tendencies (using a trapezoidal forward in time scheme)
   integer, intent(in) :: kt ! iteration index for RK (zero based - 0,1 or 2 for RK3)
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   complex(r_kind), intent(in), dimension(ndimspec,nlevs) :: &
   divspec_prev,virtempspec_prev
   complex(r_kind), intent(in), dimension(ndimspec) :: lnpsspec_prev
   real(r_double),intent(in) ::  dtx ! time step for (kt+1)'th iteration of RK
   complex(r_kind), dimension(ndimspec,nlevs) :: &
   divspec_new,virtempspec_new,espec,fspec
   complex(r_kind), dimension(ndimspec) :: lnpsspec_new,gspec
   complex(r_kind), dimension(nlevs) :: rhs
   integer n
! remove 0.5*linear terms from computed tendencies.
!$omp parallel do private(n)
   do n=1,ndimspec
      ddivspecdt(n,:) = ddivspecdt(n,:) + 0.5*lap(n)*& 
      (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      dvirtempspecdt(n,:) = dvirtempspecdt(n,:) + 0.5*matmul(bmhyb,divspec(n,:))
      dlnpsspecdt(n) = dlnpsspecdt(n) + 0.5*sum(svhyb(:)*divspec(n,:))
   enddo
!$omp end parallel do 
! solve for updated divergence.
! back subsitution to get updated virt temp, lnps.
   espec = divspec_prev + dtx*ddivspecdt
   fspec = virtempspec_prev + dtx*dvirtempspecdt
   gspec = lnpsspec_prev + dtx*dlnpsspecdt
!$omp parallel do private(n,rhs)
   do n=1,ndimspec
      rhs = espec(n,:) - 0.5*lap(n)*dtx*&
      (matmul(amhyb,fspec(n,:)) + tor_hyb(:)*gspec(n))
      divspec_new(n,:) = matmul(d_hyb_m(:,:,degree(n)+1,kt+1),rhs)
      virtempspec_new(n,:) = fspec(n,:) - 0.5*dtx*&
      matmul(bmhyb,divspec_new(n,:))
      lnpsspec_new(n) = gspec(n) - 0.5*dtx*sum(svhyb(:)*divspec_new(n,:))
   enddo
!$omp end parallel do 
! create new tendencies, including semi-implicit contribution.
   ddivspecdt = (divspec_new - divspec_prev)/dtx
   dvirtempspecdt = (virtempspec_new - virtempspec_prev)/dtx
   dlnpsspecdt = (lnpsspec_new - lnpsspec_prev)/dtx
   return
 end subroutine semimpadj

 subroutine getomega(ug,vg,divg,dlnpsdx,dlnpsdy,dlnpsdt,dlnpdtg,etadot,&
! pass in work storage so it can be re-used, saving memory.
  workb,workc,cg,cb,db)
    ! compute omega, etadot, tendency of lnps 
    ! all input and output arrays oriented bottom to top (k=1 is near ground)
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: ug,vg,divg      
    real(r_kind), intent(in), dimension(nlons,nlats) :: dlnpsdx,dlnpsdy
    ! omega (pressure vertical velocity divided by pressure) on model layers.
    real(r_kind), intent(out), dimension(nlons,nlats,nlevs) :: dlnpdtg
    ! etadot (vertical velocity in hybrid coords) on layer interfaces.
    real(r_kind), intent(out), dimension(nlons,nlats,nlevs+1) :: etadot
    real(r_kind), intent(inout), dimension(nlons,nlats) :: dlnpsdt
! work space:
    real(r_kind), intent(inout), dimension(nlons,nlats,nlevs) :: &
    workb,workc,cg,cb,db
! local scalars 
    integer k

!$omp parallel do private(k)
    do k=1,nlevs
     cg(:,:,k)=ug(:,:,nlevs+1-k)*dlnpsdx(:,:)+vg(:,:,nlevs+1-k)*dlnpsdy(:,:)
    enddo
!$omp end parallel do 
 
    db(:,:,1)=divg(:,:,nlevs)*dpk(:,:,1)
    cb(:,:,1)=cg(:,:,1)*dbk(1)
 
    do k=1,nlevs-1
     db(:,:,k+1)=db(:,:,k)+divg(:,:,nlevs-k)*dpk(:,:,k+1)
     cb(:,:,k+1)=cb(:,:,k)+cg(:,:,k+1)*dbk(k+1)
    enddo

    dlnpsdt(:,:) = -db(:,:,nlevs)/psg(:,:)-cb(:,:,nlevs)
    etadot(:,:,1) = 0.; etadot(:,:,nlevs+1) = 0.
!$omp parallel do private(k)
    do k=1,nlevs-1
       etadot(:,:,k+1)=-psg(:,:)*(bk(k+1)*dlnpsdt(:,:)+cb(:,:,k)) - db(:,:,k) 
    enddo
!$omp end parallel do 
 
    workb(:,:,1)=alfa(:,:,1)* &
                ( divg(:,:,nlevs)*dpk(:,:,1)+psg(:,:)*cb(:,:,1)*dbk(1) )
 
!$omp parallel do private(k)
    do k=2,nlevs
      workb(:,:,k)=rlnp(:,:,k)*( db(:,:,k-1)+psg(:,:)*cb(:,:,k-1) )+&
      alfa(:,:,k)*( divg(:,:,nlevs+1-k)*dpk(:,:,k)+psg(:,:)*cg(:,:,k)*dbk(k) )
    enddo
!$omp end parallel do 
 
    workc(:,:,1)=psg(:,:)*cg(:,:,1)*dbk(1)

!$omp parallel
!$omp do private(k)
    do k=2,nlevs
      workc(:,:,k)=psg(:,:)*cg(:,:,k)*( dbk(k)+ck(k)*rlnp(:,:,k)/dpk(:,:,k) )
    enddo
!$omp end do
!$omp do private(k)
    do k=1,nlevs
      dlnpdtg(:,:,nlevs+1-k)=(workc(:,:,k)-workb(:,:,k))/dpk(:,:,k)
    enddo
!$omp end do
!$omp end parallel 
 
    return
 end subroutine getomega 

 subroutine getpresgrad(virtempg,dvirtempdx,dvirtempdy,dphisdx,dphisdy,dlnpsdx,dlnpsdy,&
                        prsgx,prsgy,&
! pass in work storage so it can be re-used, saving memory.
                        cofa,cofb,px3u,px3v)
    ! compute pressure gradient terms
    ! all input and output arrays oriented bottom to top (k=1 is near ground)
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: &
    virtempg,dvirtempdx,dvirtempdy
    real(r_kind), intent(in), dimension(nlons,nlats) :: &
    dlnpsdx,dlnpsdy,dphisdx,dphisdy
    ! pressure gradient terms
    real(r_kind), intent(out), dimension(nlons,nlats,nlevs) :: &
    prsgx,prsgy
! work storage
    real(r_kind), dimension(nlons,nlats,nlevs), intent(inout) :: &
    cofa,cofb,px3u,px3v
    integer k

    cofb(:,:,1)=-(1./dpk(:,:,1))*(                 alfa(:,:,1)*dbk(1))
 
!$omp parallel
!$omp do private(k)
    do k=2,nlevs
       cofb(:,:,k)=-(1./dpk(:,:,k))*(bk(k)*rlnp(:,:,k)+alfa(:,:,k)*dbk(k))
    enddo
!$omp end do 
!$omp do private(k)
    do k=1,nlevs
        prsgx(:,:,nlevs-k+1)=cofb(:,:,k)*rd*virtempg(:,:,nlevs+1-k)*psg(:,:)*dlnpsdx(:,:)
        prsgy(:,:,nlevs-k+1)=cofb(:,:,k)*rd*virtempg(:,:,nlevs+1-k)*psg(:,:)*dlnpsdy(:,:)
    enddo
!$omp end do 
!$omp do private(k)
    do k=1,nlevs
       cofa(:,:,k)=-(1./dpk(:,:,k))*( &
        bk(k+1)*pk(:,:,k)/pk(:,:,k+1) - bk(k) &
       +rlnp(:,:,k)*( bk(k)-pk(:,:,k)*dbk(k)/dpk(:,:,k) )  )
    enddo
!$omp end do 
!$omp end parallel

    cofb(:,:,nlevs)=0.
    cofb(:,:,nlevs-1)= &
      -rd*( bk(nlevs+1)/pk(:,:,nlevs+1)-bk(nlevs)/pk(:,:,nlevs) )*virtempg(:,:,1)
 
 
    do k=2,nlevs-1
       cofb(:,:,nlevs-k)=cofb(:,:,nlevs+1-k) &
       -rd*(bk(nlevs+2-k)/pk(:,:,nlevs+2-k)-bk(nlevs+1-k)/pk(:,:,nlevs+1-k))*&
                                                        virtempg(:,:,k)
    enddo
 
    px3u(:,:,nlevs)=0. 
    px3v(:,:,nlevs)=0. 
    px3u(:,:,nlevs-1)=-rd*rlnp(:,:,nlevs)*dvirtempdx(:,:,1)
    px3v(:,:,nlevs-1)=-rd*rlnp(:,:,nlevs)*dvirtempdy(:,:,1)
    do k=2,nlevs-1
       px3u(:,:,nlevs-k)=px3u(:,:,nlevs+1-k)-rd*rlnp(:,:,nlevs+1-k)*dvirtempdx(:,:,k)
       px3v(:,:,nlevs-k)=px3v(:,:,nlevs+1-k)-rd*rlnp(:,:,nlevs+1-k)*dvirtempdy(:,:,k)
    enddo
 
!$omp parallel do private(k)
    do k=1,nlevs
       prsgx(:,:,nlevs-k+1)=prsgx(:,:,nlevs-k+1)-dphisdx(:,:)& !px1u
       +cofb(:,:,k)*psg(:,:)*dlnpsdx(:,:)& !px2u
       +px3u(:,:,k)&
       -rd*alfa(:,:,k)*dvirtempdx(:,:,nlevs+1-k)& !px4u
       -cofa(:,:,k)*rd*virtempg(:,:,nlevs+1-k)*psg(:,:)*dlnpsdx(:,:) !px5u
       prsgy(:,:,nlevs-k+1)=prsgy(:,:,nlevs-k+1)-dphisdy(:,:)& !px1v
       +cofb(:,:,k)*psg(:,:)*dlnpsdy(:,:)& !px2v
       +px3v(:,:,k)&
       -rd*alfa(:,:,k)*dvirtempdy(:,:,nlevs+1-k)& !px4v
       -cofa(:,:,k)*rd*virtempg(:,:,nlevs+1-k)*psg(:,:)*dlnpsdy(:,:) !px5v
    enddo
!$omp end parallel do 

 end subroutine getpresgrad

 subroutine getvadv(datag,etadot,vadv)
   ! compute vertical advection of datag, using etadot - result in vadv
   ! datag, vadv bottom to top, etadot top to bottom.
   real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: datag
   real(r_kind), intent(in), dimension(nlons,nlats,nlevs+1) :: etadot
   real(r_kind), intent(out), dimension(nlons,nlats,nlevs) :: vadv
   integer k

   vadv(:,:,nlevs)= &
   (0.5/dpk(:,:,1))*etadot(:,:,2)*( datag(:,:,nlevs-1)-datag(:,:,nlevs))

   vadv(:,:,1)= &
   (0.5/dpk(:,:,nlevs))*etadot(:,:,nlevs)*( datag(:,:,1)-datag(:,:,2) )
 
!$omp parallel do private(k)
   do k=2,nlevs-1
      vadv(:,:,nlevs+1-k)= &
   (0.5/dpk(:,:,k))*( etadot(:,:,k+1)*( datag(:,:,nlevs  -k)-datag(:,:,nlevs+1-k) )+&
                      etadot(:,:,k  )*( datag(:,:,nlevs+1-k)-datag(:,:,nlevs+2-k) ) )
   enddo
!$omp end parallel do 

   return
 end subroutine getvadv

 subroutine getvadv_tracers(datag,etadot,vadv)
   ! compute vertical advection of datag, using etadot - result in vadv
   ! datag, vadv bottom to top, etadot top to bottom.
   ! uses positive definite scheme of Thuburn (1993, DOI:
   ! 10.1002/qj.49711951107)
   real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: datag
   real(r_kind), intent(in), dimension(nlons,nlats,nlevs+1) :: etadot
   real(r_kind), intent(out), dimension(nlons,nlats,nlevs) :: vadv
   real(r_kind), dimension(:,:,:), allocatable :: datag_half, datag_d
   integer i,j,k
   real(r_kind) rrkp,rrk1m,phkp,phkp1m

   allocate(datag_half(nlons,nlats,0:nlevs))
   allocate(datag_d(nlons,nlats,0:nlevs))

!$omp parallel do private(k)
   do k=1,nlevs-1
      datag_half(:,:,k) = 0.5*(datag(:,:,nlevs-k)+datag(:,:,nlevs+1-k))
   enddo
!$omp end parallel do 
   datag_half(:,:,0) = datag(:,:,nlevs)
   datag_half(:,:,nlevs) = datag(:,:,1)
!$omp parallel do private(k)
   do k=1,nlevs-1
      datag_d(:,:,k) = datag(:,:,nlevs-k) - datag(:,:,nlevs+1-k)
   enddo
!$omp end parallel do 

!$omp parallel do private(i,j)
   do j=1,nlats
   do i=1,nlons
      if (datag(i,j,nlevs) >= 0.) then
         datag_d(i,j,0) = datag(i,j,nlevs) - &
         max(0.,2.*datag(i,j,nlevs)-datag(i,j,nlevs-1))
      else
         datag_d(i,j,0) = datag(i,j,nlevs) - &
         min(0.,2.*datag(i,j,nlevs)-datag(i,j,nlevs-1))
      end if
      if (datag(i,j,1) >= 0) then
         datag_d(i,j,nlevs) = max(0.,2.*datag(i,j,1)-datag(i,j,2)) - datag(i,j,1)
      else
         datag_d(i,j,nlevs) = min(0.,2.*datag(i,j,1)-datag(i,j,2)) - datag(i,j,1)
      endif
   enddo
   enddo
!$omp end parallel do 

!$omp parallel do private(i,j,k,rrkp,phkp,rrk1m,phkp1m)
   do k=1,nlevs-1
      do j=1,nlats
      do i=1,nlons
         if(etadot(i,j,k+1) > 0.) then  !etadot is from top to bottom
            rrkp = 0.
            if (datag_d(i,j,k) .ne. 0.) rrkp = datag_d(i,j,k-1)/datag_d(i,j,k)
            phkp = (rrkp+abs(rrkp))/(1.+abs(rrkp))
            datag_half(i,j,k) = datag(i,j,nlevs+1-k) + &
                                phkp*(datag_half(i,j,k)-datag(i,j,nlevs+1-k))
         else
            rrk1m = 0.
            if (datag_d(i,j,k) .ne. 0.) rrk1m = datag_d(i,j,k+1)/datag_d(i,j,k)
            phkp1m = (rrk1m+abs(rrk1m))/(1.+abs(rrk1m))
            datag_half(i,j,k) = datag(i,j,nlevs-k) + &
                                phkp1m*(datag_half(i,j,k)-datag(i,j,nlevs-k))
         endif
      enddo
      enddo
   enddo
!$omp end parallel do 

!$omp parallel do private(k)
   do k=1,nlevs
      vadv(:,:,nlevs+1-k) = (1./dpk(:,:,k))*(&
      (datag_half(:,:,k)*etadot(:,:,k+1) - datag_half(:,:,k-1)*etadot(:,:,k))+&
      (datag(:,:,nlevs+1-k)*(etadot(:,:,k)-etadot(:,:,k+1))))
   enddo
!$omp end parallel do 
      
   deallocate(datag_half, datag_d)

   return
 end subroutine getvadv_tracers

end module dyn_run
