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

 use params, only: nlevs,ntrunc,nlons,nlats,ndimspec,dt,ntrac,pdryini,&
   dcmip,dry,vcamp,svc
 use kinds, only: r_kind,r_double
 use shtns, only: degree,order,&
 lap,invlap,lats,grdtospec,spectogrd,getuv,getvrtdivspec,getgrad,areawts
 use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
 tracerspec, topospec
 use grid_data, only: ug,vg,vrtg,virtempg,divg,tracerg,dlnpdtg,etadot,lnpsg,&
 phis,dphisdx,dphisdy,dlnpsdt,keg
 use pressure_data, only:  ak, bk, ck, dbk, dpk, rlnp, pk, alfa, dpk, psg,&
 calc_pressdata
 use stoch_data, only:  grd_svc, vfact_svc
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
               eps => con_eps, omega => con_omega, cvap => con_cvap, &
               grav => con_g, fv => con_fvirt, kappa => con_rocp
 use semimp_data, only: amhyb,bmhyb,svhyb,d_hyb_m,tor_hyb,tref,pkref,alfaref,&
                        getdtfact

 implicit none
 private

 public :: getdyntend, getomega, getpresgrad, getvadv, semimpadj

 contains

 subroutine getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,&
                       dtracerspecdt,dlnpsspecdt,kstep,just_do_inverse_transform)
   ! compute dynamics tendencies (not including hyper-diffusion.and linear drag,
   ! which are treated implicitly)
   ! based on hybrid sigma-pressure dynamical core described in
   ! http://www.emc.ncep.noaa.gov/officenotes/newernotes/on462.pdf
   ! The only difference is that here we use a forward in time
   ! runge-kutta scheme instead of robert-assellin time filtered leap-frog.
   logical, optional, intent(in) :: just_do_inverse_transform
   logical :: early_return = .false. ! if true. spectral-> grid only
   complex(r_kind), intent(out), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
   dtracerspecdt
   complex(r_kind), intent(out), dimension(ndimspec) :: dlnpsspecdt
! local variables.
   complex(r_kind), dimension(:,:),allocatable :: workspec
   real(r_kind), dimension(ndimspec) :: smoothfact
   real(r_kind), dimension(:,:), allocatable :: dlnpsdx,dlnpsdy
   real(r_kind), dimension(:,:,:), allocatable :: &
   prsgx,prsgy,vadvu,vadvv,vadvt,vadvq,dvirtempdx,dvirtempdy
   integer, intent(in) :: kstep ! runge-kutta step (0,1,2)
   integer k,nt,ntrac_use
   logical :: profile = .false. ! print out timing stats
   real(r_kind) pmoist,pdry,epstiny,rnn1,rnn0
   real(8) t1,t2,t0
   integer(8) count, count_rate, count_max

   epstiny = tiny(epstiny)

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
   if (vcamp > epstiny .or. svc > epstiny) then
      ! smoothing factors used to compute vorticity confinement
      ! anti-diffusion of vorticity.
      rnn0 = 0.5*ntrunc*(0.5*ntrunc+1)
      do k=1,ndimspec
         rnn1 = degree(k)*(degree(k)+1)
         smoothfact(k)=exp(-(rnn1/rnn0))
      enddo
   endif

   ! compute u,v,virt temp, vorticity, divergence, ln(ps)
   ! and specific humidity on grid from spectral coefficients.
   call system_clock(count, count_rate, count_max)
   t1 = count*1.d0/count_rate
   t0 = t1
!$omp parallel do private(k,nt) schedule(dynamic)
   do k=1,nlevs
      call getuv(vrtspec(:,k),divspec(:,k),ug(:,:,k),vg(:,:,k),rerth)
      call spectogrd(vrtspec(:,k),vrtg(:,:,k))
      call spectogrd(divspec(:,k),divg(:,:,k))
      call spectogrd(virtempspec(:,k),virtempg(:,:,k))
      ! gradient of virtual temperature on grid.
      call getgrad(virtempspec(:,k),dvirtempdx(:,:,k),dvirtempdy(:,:,k),rerth)
      ! specific humidity, other tracers on grid.
      do nt=1,ntrac
         call spectogrd(tracerspec(:,k,nt),tracerg(:,:,k,nt))
      enddo
   enddo
!$omp end parallel do 
   call system_clock(count, count_rate, count_max)
   t2 = count*1.d0/count_rate
   if (profile) print *,'1 time=',t2-t1

   ! lnps on grid.
   call spectogrd(lnpsspec, lnpsg)

   ! compute pressure on interfaces and model layers using lnpsg
   ! i.e. calculate alpha,delp,rlnp,dpk,psg,pk from ak,bk,lnpsg
   ! results stored in module pressure_data
   call calc_pressdata(lnpsg) 

! compute global mean dry ps (only first step in RK scheme).
   if (kstep .eq. 0) then
      pmoist = 0.
      if (.not. dry) then
         ntrac_use = ntrac
         if (dcmip > 0) ntrac_use = 1
         do k=1,nlevs
            do nt=1,ntrac_use
               pmoist = pmoist + sum(areawts*dpk(:,:,nlevs-k+1)*tracerg(:,:,k,nt))
            enddo
         enddo
      endif
      pdry = sum(areawts*psg) - pmoist
      if (pdryini .eq. 0) then
         pdryini = pdry
         print *,'pdryini = ',pdryini
      else
         print *,'pdry = ',pdry
      endif
   endif

   ! gradient of surface pressure.
   ! gradient of surface geopotential precomputed in dyn_init,
   ! saved in module grid_data (dphisdx,dphisdy).
   call getgrad(lnpsspec, dlnpsdx, dlnpsdy, rerth)

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

   ! get tendency of lnps on grid.
   call grdtospec(dlnpsdt,dlnpsspecdt)
   ! get pressure gradient terms (prsgx,prsgy)
   call getpresgrad(virtempg,dvirtempdx,dvirtempdy,dphisdx,dphisdy,dlnpsdx,dlnpsdy,&
                    prsgx,prsgy,&
                    vadvu,vadvv,vadvt,vadvq) ! work storage

   ! get vertical advection terms  (input etadot is top to bottom)
   call getvadv(ug,etadot,vadvu)
   call getvadv(vg,etadot,vadvv)
   call getvadv(virtempg,etadot,vadvt)

   ! temporary storage of planetary vorticity
   dlnpsdx = 2.*omega*sin(lats) 

   ! compute energy conversion term, store in vadvq.
   if (ntrac > 0) then
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
   if (profile) print *,'2 time=',t1-t2

   ! compute tendencies of virt temp, vort and div in spectral space
!$omp parallel do private(k) schedule(dynamic)
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
      if (vcamp > epstiny .or. svc > epstiny) then
         ! add simplified vorticity confinement (anti-diffusion of vorticity)
         ! operates on smoothed vorticity field.
         ! multiply anti-diffusion by stochastic random pattern if svc > 0
         if (svc > epstiny) then
            workspec(:,k) = (lap/rerth**2)*smoothfact*vrtspec(:,k)
            call spectogrd(workspec(:,k), prsgx(:,:,k))
            prsgx(:,:,k) = (vcamp+vfact_svc(k)*grd_svc)*prsgx(:,:,k)
            call grdtospec(prsgx(:,:,k),workspec(:,k))
            dvrtspecdt(:,k) = dvrtspecdt(:,k) - workspec(:,k)
         else ! if not stochastic, no transforms needed (VC is cost-free)
            dvrtspecdt(:,k) = dvrtspecdt(:,k) - vcamp*(lap/rerth**2)*smoothfact*vrtspec(:,k)
         endif
      endif
      ! add laplacian(KE) term to div tendency
      keg(:,:,k) = 0.5*(ug(:,:,k)**2+vg(:,:,k)**2)
      call grdtospec(keg(:,:,k),workspec(:,k))
      ddivspecdt(:,k) = ddivspecdt(:,k) - &
      (lap(:)/rerth**2)*workspec(:,k)
   enddo
!$omp end parallel do 
   call system_clock(count, count_rate, count_max)
   t2 = count*1.d0/count_rate
   if (profile) print *,'3 time=',t2-t1

   ! compute tendency of tracers (including specific humidity) in spectral space.
   do nt=1,ntrac
   ! use positive-definite vertical advection.
   call getvadv_tracers(tracerg(:,:,:,nt),etadot,vadvq)
!$omp parallel do private(k) schedule(dynamic)
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
   if (profile) print *,'4 time=',t1-t2  
   !print *,'getdyntend time=',t1-t0

   ! deallocate storage.
   deallocate(vadvq,workspec,dvirtempdx,dvirtempdy)
   deallocate(prsgx,prsgy,vadvu,vadvv,vadvt)
   deallocate(dlnpsdx,dlnpsdy)

   return
 end subroutine getdyntend

 subroutine semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                      divspec_orig,virtempspec_orig,lnpsspec_orig,kt,dtx)
! semi-implicit adjustment of tendencies (using a trapezoidal forward in time scheme)
   integer, intent(in) :: kt ! RK stage (zero based - 0,1 or 2 for 3-stage RK)
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   complex(r_kind), intent(in), dimension(ndimspec,nlevs) :: &
   divspec_orig,virtempspec_orig
   complex(r_kind), intent(in), dimension(ndimspec) :: lnpsspec_orig
   real(r_double),intent(in) ::  dtx ! time step for (kt+1)'th RK stage.
   complex(r_kind) :: espec(nlevs),fspec(nlevs),gspec,rhs(nlevs),&
    divav(nlevs),vtempav(nlevs),lnpsav
   integer n
   real(r_kind) dtfact_new,dtfact_orig,dtfact_current

   ! dtfact_new is coeff for next stage, dtfact_current is coeff for current
   ! stage, dtfact_orig is coeff for value at beginning of RK cycle.
   call getdtfact(kt,dtfact_orig,dtfact_current,dtfact_new)

! solve for updated divergence.
! back subsitution to get updated virt temp, lnps.

!$omp parallel do private(n,rhs,espec,fspec,gspec,divav,lnpsav,vtempav)
   do n=1,ndimspec
! first remove linear terms from computed tendencies.
      ddivspecdt(n,:) = ddivspecdt(n,:) + lap(n)*& 
      (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      dvirtempspecdt(n,:) = dvirtempspecdt(n,:) + matmul(bmhyb,divspec(n,:))
      dlnpsspecdt(n) = dlnpsspecdt(n) + sum(svhyb(:)*divspec(n,:))
! solve for updated divergence
      vtempav  = dtfact_current*virtempspec(n,:)+dtfact_orig*virtempspec_orig(n,:)
      divav    = dtfact_current*divspec(n,:)    +dtfact_orig*divspec_orig(n,:)
      lnpsav   = dtfact_current*lnpsspec(n)     +dtfact_orig*lnpsspec_orig(n)
      espec(:) = divspec_orig(n,:) + dtx*ddivspecdt(n,:) - & 
                 dtx*lap(n)*(matmul(amhyb,vtempav(:)) + &
                                    tor_hyb(:)*lnpsav)
      fspec(:) = virtempspec_orig(n,:) + dtx*dvirtempspecdt(n,:) - &
                 dtx*matmul(bmhyb, divav(:))
      gspec    = lnpsspec_orig(n) + dtx*dlnpsspecdt(n) - &
                 dtx*sum(svhyb(:)*divav(:))
      rhs(:)   = espec(:) - dtfact_new*lap(n)*dtx*&
                 (matmul(amhyb,fspec(:)) + tor_hyb(:)*gspec)
! back substitution to get updated virt temp, lnps.
! create new tendencies, including semi-implicit contribution.
      ! espec holds updated divergence.
      espec(:) = matmul(d_hyb_m(:,:,degree(n)+1,kt+1),rhs(:))
      ddivspecdt(n,:) = (espec(:)-divspec_orig(n,:))/dtx
      dvirtempspecdt(n,:) = (fspec(:) - dtfact_new*dtx*matmul(bmhyb,espec(:))-&
                             virtempspec_orig(n,:))/dtx
      dlnpsspecdt(n) = (gspec - dtfact_new*dtx*sum(svhyb(:)*espec(:))-&
                        lnpsspec_orig(n))/dtx
   enddo
!$omp end parallel do 

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
   ! local variables 
   real(r_kind), dimension(:,:,:), allocatable :: datag_half, datag_d
   integer i,j,k
   real(r_kind) epstiny,phi
   !real(r_kind) phi(nlons,nlats)

   allocate(datag_half(nlons,nlats,0:nlevs))
   allocate(datag_d(nlons,nlats,0:nlevs))

!$omp parallel do private(k)
   do k=1,nlevs-1
      datag_half(:,:,k) = 0.5*(datag(:,:,nlevs-k)+datag(:,:,nlevs+1-k))
      datag_d(:,:,k) = datag(:,:,nlevs-k) - datag(:,:,nlevs+1-k)
   enddo
!$omp end parallel do 

   datag_half(:,:,0) = datag(:,:,nlevs)
   datag_half(:,:,nlevs) = datag(:,:,1)

!$omp parallel do private(i)
   do i=1,nlons
      where (datag(i,:,nlevs) >= 0)
         datag_d(i,:,0) = datag(i,:,nlevs) - &
         max(0.,2.*datag(i,:,nlevs)-datag(i,:,nlevs-1))
      elsewhere
         datag_d(i,:,0) = datag(i,:,nlevs) - &
         min(0.,2.*datag(i,:,nlevs)-datag(i,:,nlevs-1))
      end where
      where (datag(i,:,1) >= 0)
         datag_d(i,:,nlevs) = max(0.,2.*datag(i,:,1)-datag(i,:,2)) - datag(i,:,1)
      elsewhere
         datag_d(i,:,nlevs) = min(0.,2.*datag(i,:,1)-datag(i,:,2)) - datag(i,:,1)
      end where
   enddo
!$omp end parallel do 

! same as above, but without where blocks inside parallel region.
!!$omp parallel do private(i,j)
!   do i=1,nlons
!      do j=1,nlats
!         if (datag(i,j,nlevs) >= 0) then
!            datag_d(i,j,0) = datag(i,j,nlevs) - &
!            max(0.,2.*datag(i,j,nlevs)-datag(i,j,nlevs-1))
!         else
!            datag_d(i,j,0) = datag(i,j,nlevs) - &
!            min(0.,2.*datag(i,j,nlevs)-datag(i,j,nlevs-1))
!         endif
!         if (datag(i,j,1) >= 0) then
!            datag_d(i,j,nlevs) = max(0.,2.*datag(i,j,1)-datag(i,j,2)) - datag(i,j,1)
!         else
!            datag_d(i,j,nlevs) = min(0.,2.*datag(i,j,1)-datag(i,j,2)) - datag(i,j,1)
!         endif
!      enddo
!   enddo
!!$omp end parallel do 

   ! to prevent NaNs in computation of van leer limiter
   epstiny = tiny(epstiny)
!$omp parallel do private(i,j,k)
   do k=1,nlevs-1
      do j=1,nlats
      do i=1,nlons
         if (abs(datag_d(i,j,k)) .lt. epstiny) then
            datag_d(i,j,k) = sign(epstiny, datag_d(i,j,k))
         endif
      enddo
      enddo
   enddo
!$omp end parallel do 

!   ! apply flux limiter.
!!note: this segfaults with intel 12.1 at T574
!!$omp parallel do private(k,phi)
!   do k=1,nlevs-1
!      where(etadot(:,:,k+1) > 0.)
!         phi = datag_d(:,:,k-1)/datag_d(:,:,k)
!         datag_half(:,:,k) = datag(:,:,nlevs+1-k) + &
!         (phi+abs(phi))/(1.+abs(phi))*(datag_half(:,:,k)-datag(:,:,nlevs+1-k))
!      elsewhere
!         phi = datag_d(:,:,k+1)/datag_d(:,:,k)
!         datag_half(:,:,k) = datag(:,:,nlevs-k) + &
!         (phi+abs(phi))/(1.+abs(phi))*(datag_half(:,:,k)-datag(:,:,nlevs-k))
!      end where
!   enddo
!!$omp end parallel do 

!    ! apply flux limiter.
! to avoid t574 segfault, get rid of where statements inside parallel region.
!$omp parallel do private(i,j,k,phi)
    do k=1,nlevs-1
       do j=1,nlats
       do i=1,nlons
          if (etadot(i,j,k+1) > 0.) then
             phi = datag_d(i,j,k-1)/datag_d(i,j,k)
             datag_half(i,j,k) = datag(i,j,nlevs+1-k) + &
             (phi+abs(phi))/(1.+abs(phi))*(datag_half(i,j,k)-datag(i,j,nlevs+1-k))
          else
             phi = datag_d(i,j,k+1)/datag_d(i,j,k)
             datag_half(i,j,k) = datag(i,j,nlevs-k) + &
             (phi+abs(phi))/(1.+abs(phi))*(datag_half(i,j,k)-datag(i,j,nlevs-k))
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
