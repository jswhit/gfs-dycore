module stoch_data

! initialize data structures and arrays for stochastic perturbation physics.

 use kinds, only: r_kind, r_single
 use params, only: dt,svc,sppt,ndimspec,nlons,nlats,ntrunc,nlevs,&
 svc_tau,svc_lscale,sppt_tau,sppt_lscale,&
 iseed_svc,iseed_sppt,iseed_shum,shum,shum_tau,shum_lscale,&
 iseed_addnoise,addnoise,addnoise_tau,addnoise_lscale,addnoise_vfilt,&
 addnoise_vrtonly
 use patterngenerator, only: random_pattern, patterngenerator_init,&
 getnoise, patterngenerator_advance, patterngenerator_destroy
 use physcons, only: pi => con_pi
 use pressure_data, only: sl

 implicit none

 private

 public :: init_stochdata,destroy_stochdata,getstochforcing
 complex(r_kind), allocatable, public, dimension(:) :: &
 spec_svc,spec_sppt,spec_shum,specps_addnoise
 complex(r_kind), allocatable, public, dimension(:,:) :: &
 specpsi_addnoise,spect_addnoise,specchi_addnoise
 real(r_kind), allocatable, public, dimension(:,:) :: &
 grd_svc,grd_sppt,grd_shum
 real(r_kind), allocatable, public, dimension(:) :: &
 vfact_shum,vfact_svc,vfact_sppt,vfact_addnoise
 real(r_kind), allocatable, dimension(:,:,:) :: agv
 type(random_pattern), public :: rpattern_svc
 type(random_pattern), public :: rpattern_sppt
 type(random_pattern), public :: rpattern_shum
 type(random_pattern), public :: rpattern_addnoise

 contains

 subroutine init_stochdata()

    real(r_kind), allocatable, dimension(:,:,:) :: agv_tmp
    real(r_single), allocatable, dimension(:,:,:) :: agvin
    real(r_single), allocatable, dimension(:,:) :: bvin,wgvin
    real(r_kind), allocatable, dimension(:) :: akin,bkin,siin,rlsigo,rlsig,&
      coef1,coef2
    integer,allocatable, dimension(:) :: lsig
    integer n,nspinup,k,nlatsin,nlevsin,i,j,l,l1,m,m1,spinup_efolds

    real(r_kind) delt,sigtop,sigbot
    delt = dt
    spinup_efolds = 5
    if (svc > tiny(svc)) then
! vorticity confinement.
       allocate(spec_svc(ndimspec))
       allocate(grd_svc(nlons,nlats))
       allocate(vfact_svc(nlevs))
       call patterngenerator_init(svc_lscale,delt,svc_tau,svc,iseed_svc,rpattern_svc,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = spinup_efolds*svc_tau/delt
       call getnoise(rpattern_svc,spec_svc)
       spec_svc = rpattern_svc%stdev*spec_svc*rpattern_svc%varspectrum
       do n=1,nspinup
          call patterngenerator_advance(spec_svc,rpattern_svc)
       enddo
       sigbot = 0.10; sigtop = 0.05
       do k=1,nlevs
          vfact_svc(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigtop) then
            vfact_svc(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_svc(k) = 0.
          endif
       enddo
    endif
    if (sppt > tiny(sppt)) then
! stochastically perturbed physics tendencies.
       allocate(spec_sppt(ndimspec))
       allocate(grd_sppt(nlons,nlats))
       allocate(vfact_sppt(nlevs))
       call patterngenerator_init(sppt_lscale,delt,sppt_tau,sppt,iseed_sppt,rpattern_sppt,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = spinup_efolds*sppt_tau/delt
       call getnoise(rpattern_sppt,spec_sppt)
       spec_sppt = rpattern_sppt%stdev*spec_sppt*rpattern_sppt%varspectrum
       do n=1,nspinup
          call patterngenerator_advance(spec_sppt,rpattern_sppt)
       enddo
       sigbot = 0.10; sigtop = 0.05
       do k=1,nlevs
          vfact_sppt(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigtop) then
            vfact_sppt(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_sppt(k) = 0.
          endif
       enddo
    endif
    if (shum > tiny(shum)) then
! perturbed boundary layer specific humidity.
       allocate(spec_shum(ndimspec))
       allocate(grd_shum(nlons,nlats))
       allocate(vfact_shum(nlevs))
       call patterngenerator_init(shum_lscale,delt,shum_tau,shum,iseed_shum,rpattern_shum,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = spinup_efolds*shum_tau/delt
       call getnoise(rpattern_shum,spec_shum)
       spec_shum = rpattern_shum%stdev*spec_shum*rpattern_shum%varspectrum
       sigbot = 0.2
       ! humidity pert decays exponentially away from sfc
       do k=1,nlevs
          vfact_shum(k) = exp((sl(k)-1.)/sigbot)
       enddo
       do n=1,nspinup
          call patterngenerator_advance(spec_shum,rpattern_shum)
       enddo
    endif
    if (addnoise > tiny(addnoise)) then
! additive noise.
       allocate(specpsi_addnoise(ndimspec,nlevs))
       allocate(vfact_addnoise(nlevs))

       call patterngenerator_init(addnoise_lscale,delt,addnoise_tau,&
                                  addnoise,iseed_addnoise,rpattern_addnoise,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = spinup_efolds*addnoise_tau/delt
       do k=1,nlevs
          call getnoise(rpattern_addnoise,specpsi_addnoise(:,k))
          specpsi_addnoise(:,k) = &
          rpattern_addnoise%stdev*specpsi_addnoise(:,k)*rpattern_addnoise%varspectrum(:)
       enddo
       sigbot = 0.10; sigtop = 0.05
       do k=1,nlevs
          vfact_addnoise(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigtop) then
            vfact_addnoise(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_addnoise(k) = 0.
          endif
       enddo
       do n=1,nspinup
       do k=1,nlevs
          call patterngenerator_advance(specpsi_addnoise(:,k),rpattern_addnoise)
       enddo
       enddo
! if not addnoise_vrtonly, read in balance operators from GSI to
! create perts for tv.
       if (.not. addnoise_vrtonly) then
!         read in balance operators from GSI berror_stats file.
          print *,'read in bal operators'
          open(7,file='berror_stats',form='unformatted',status='old',convert='big_endian')
          read(7) nlevsin,nlatsin
          if (nlatsin .ne. nlats+2) then
             print *,'incorrect # of lats in berror_stats',nlatsin
             stop
          endif
          allocate ( agvin(nlats+2,nlevsin,nlevsin) )
          allocate ( bvin(nlats+2,nlevsin),wgvin(nlats+2,nlevsin) )
          allocate ( agv_tmp(nlats,nlevsin,nlevsin) )
          allocate ( agv(nlats,nlevs,nlevs) )
          !   agv,wgv,bv - balance regression matrix for t,ps,div
          read(7) agvin,bvin,wgvin
          close(7)
          do k=1,nlevsin
            do n=1,nlevsin
               agv_tmp(:,k,n) = agvin(2:nlats+1,k,n)
            enddo
          enddo
!   compute vertical(pressure) interpolation index and weight
          allocate(akin(nlevsin+1),bkin(nlevsin+1))
          allocate(siin(nlevsin+1),rlsigo(nlevsin),rlsig(nlevs))
          print *,'interpolate balance operators'
          open(8,form='formatted',file='berror_levs')
          read(8,*) k,n
          if (n .ne. nlevsin) then
            print *,'wrong number of levs in berror_levs'
            stop
          endif
          do k=1,nlevsin+1
             read(8,*) akin(k),bkin(k)
             siin(k) = akin(k)/101300.0+bkin(k)  
          enddo
          close(8)
          do k=1,nlevsin
             rlsigo(k) = -log(0.5*(siin(k)+siin(k+1)))
             !print *,k,rlsigo(k)
          enddo
          do k=1,nlevs
             rlsig(k)=-log(sl(k))
             !print *,k,rlsig(k)
          enddo

          allocate(lsig(nlevs),coef1(nlevs),coef2(nlevs))
          do k=1,nlevs
             if(rlsig(k)<=rlsigo(1))then
                m=1
                m1=2
                lsig(k)=1
                coef1(k)=1.
                coef2(k)=0
             else if(rlsig(k)>=rlsigo(nlevsin))then
                m=nlevsin-1
                m1=nlevsin
                lsig(k)=nlevsin-1
                coef1(k)=0.
                coef2(k)=1.
             else
                m_loop: do m=1,nlevsin-1
                   m1=m+1
                   if((rlsig(k)>=rlsigo(m))   .and.  &
                        (rlsig(k)<rlsigo(m1))     )then
                      lsig(k)=m
                      exit m_loop
                   end if
                enddo m_loop
                coef1(k)=(rlsigo(m1)-rlsig(k))/(rlsigo(m1)-rlsigo(m))
                coef2(k)=1.-coef1(k)
                if(lsig(k)==nlevsin)then
                   lsig(k)=nlevsin-1
                   coef2(k)=1.
                   coef1(k)=0.
                endif
             endif
             !print *,k,lsig(k),coef1(k),coef2(k)
          end do

!   Load agv, interpolate to model levels.
          do k=1,nlevs
             m=lsig(k)
             m1=m+1
             do j=1,nlevs
                l=lsig(j)
                l1=l+1
                do i=1,nlats
                   agv(i,j,k)=(agv_tmp(i,l,m)*coef1(j)+agv_tmp(i,l1,m)*coef2(j))*coef1(k) &
                             +(agv_tmp(i,l,m1)*coef1(j)+agv_tmp(i,l1,m1)*coef2(j))*coef2(k)
                enddo
             enddo
          enddo

          deallocate(akin,bkin,rlsigo,rlsig,siin,lsig,coef1,coef2)
          deallocate(agv_tmp,agvin,wgvin,bvin)
       endif ! .not. addnoise_vrtonly
    endif ! addnoise > 0
 end subroutine init_stochdata

 subroutine destroy_stochdata()
    if (allocated(spec_svc)) then
        deallocate(spec_svc,vfact_svc,grd_svc)
        call patterngenerator_destroy(rpattern_svc)
    endif
    if (allocated(spec_sppt)) then
        deallocate(spec_sppt,vfact_sppt,grd_sppt)
        call patterngenerator_destroy(rpattern_sppt)
    endif
    if (allocated(spec_shum)) then
        deallocate(spec_shum,vfact_shum,grd_shum)
        call patterngenerator_destroy(rpattern_shum)
    endif
    if (allocated(specpsi_addnoise)) then
        deallocate(specpsi_addnoise,vfact_addnoise)
        call patterngenerator_destroy(rpattern_addnoise)
    endif
 end subroutine destroy_stochdata
 
 subroutine getstochforcing(specpsi_addnoise,psiforcing,tvforcing)
     use grid_data, only: ug,vg
     use params, only: addnoise_vfilt, addnoise_dissfact
     use physcons, only: rerth => con_rerth, pi => con_pi
     use shtns, only: degree, lap, getuv, grdtospec, spectogrd, lats, lons
     use spectral_data, only: vrtspec, divspec, disspec
     real(r_kind), dimension(nlons,nlats,nlevs) :: udiffg,vdiffg,workg
     real(r_kind), dimension(ndimspec) :: smoothfact
     complex(r_kind), dimension(ndimspec,nlevs) :: workspec,workspec2
     complex(r_kind), intent(out), dimension(ndimspec,nlevs) :: &
         psiforcing,tvforcing
     complex(r_kind), intent(in), dimension(ndimspec,nlevs) :: specpsi_addnoise
     integer n,k,i,j,l
     real(r_kind) rnn0,rnn1,epstiny
     epstiny = tiny(rnn0)
     ! apply successive applications of 1-2-1 filter in vertical to introduce vertical correlations.
     if (addnoise_vfilt > 0) then
        do n=1,addnoise_vfilt
           do k=2,nlevs-1
              workspec(:,k) = specpsi_addnoise(:,k+1)+2.*specpsi_addnoise(:,k)+&
                              specpsi_addnoise(:,k-1)
           enddo
           workspec(:,1) = (1.+1./3.)*specpsi_addnoise(:,2)+2.*(1.+1./3.)*specpsi_addnoise(:,1)
           workspec(:,nlevs) = (1.+1./3.)*specpsi_addnoise(:,nlevs-1)+2.*(1.+1./3.)*specpsi_addnoise(:,nlevs)
           psiforcing = 0.25*workspec
        enddo
     else
        psiforcing = specpsi_addnoise
     end if
     if (addnoise_dissfact) then
        ! gaussian smoothing parameter.
        ! filter length scale three times the scale of the additive noise
        ! (using weaver and courtier DOI: 10.1002/qj.49712757518)
        rnn1 = 2.*rerth/(3.*addnoise_lscale)
        print *,'rnn1 = ',rnn1
        rnn0 = rnn1*(rnn1+1.)
        do k=1,ndimspec
           rnn1 = degree(k)*(degree(k)+1)
           smoothfact(k)=exp(-(rnn1/rnn0))
        enddo
!$omp parallel do private(k)
        do k=1,nlevs
           workspec(:,k) = -disspec(:)*vrtspec(:,k)
           workspec2(:,k)= -disspec(:)*divspec(:,k)
           call getuv(workspec(:,k),workspec2(:,k),udiffg(:,:,k),vdiffg(:,:,k),rerth)
           udiffg(:,:,k) = abs(ug(:,:,k)*udiffg(:,:,k)+vg(:,:,k)*vdiffg(:,:,k))
           ! clip neg values, take square root.
           where (udiffg(:,:,k) < epstiny) udiffg(:,:,k) = epstiny
           udiffg(:,:,k) = sqrt(udiffg(:,:,k))
           ! smoothing in spectral space.
           call grdtospec(udiffg(:,:,k),workspec(:,k))
           workspec(:,k) = smoothfact*workspec(:,k)
           call spectogrd(workspec(:,k),vdiffg(:,:,k))
           !vdiffg(:,:,k) = sqrt(udiffg(:,:,k))
           !print *,k,minval(vdiffg(:,:,k)),maxval(vdiffg(:,:,k))
           call spectogrd(psiforcing(:,k),workg(:,:,k))
           ! modulate additive noise by ke dissipation rate
           udiffg(:,:,k) = vfact_addnoise(k)*workg(:,:,k)*vdiffg(:,:,k)
           !udiffg(:,:,k) =& ! for testing
           !200.*(sin(2.*lats)**6)*(sin(0.5*lons-0.5*pi)**8)*exp(-((sl(k)-0.5)/0.2)**2)
        enddo
!$omp end parallel do
     endif
     print *,'min/max modulated psi forcing',minval(udiffg),maxval(udiffg)
     print *,'min/max orig psi forcing',minval(workg),maxval(workg)
     print *,'min/max smoothed ke diss',minval(vdiffg),maxval(vdiffg)
     !open(7,form='unformatted',file='addnoise.dat',access='direct',recl=nlons*nlats*nlevs)
     !write(7,rec=1) udiffg
     !write(7,rec=2) vdiffg
     !write(7,rec=3) workg
! create balanced temp forcing (if addnoise_vrtonly is false).
     vdiffg = 0.
!$omp parallel do private(i,j,l,k)
     do k=1,nlevs
        if (.not. addnoise_vrtonly) then
           do l=1,nlevs
              do j=1,nlons
                 do i=1,nlats
                    vdiffg(j,i,k)=vdiffg(j,i,k)+agv(i,k,l)*udiffg(j,i,l)
                 end do
              end do
           enddo
           call grdtospec(vdiffg(:,:,k),tvforcing(:,k))
        else
           tvforcing(:,k) = 0.
        endif
        ! psi forcing (in udiffg) in spectral space.
        call grdtospec(udiffg(:,:,k),psiforcing(:,k))
     end do
!$omp end parallel do
     !write(7,rec=4) vdiffg
     !close(7)
     !stop
 end subroutine getstochforcing

end module stoch_data
