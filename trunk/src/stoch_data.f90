module stoch_data

 use kinds, only: r_kind
 use params, only: dt,svc,sppt,ndimspec,nlons,nlats,ntrunc,nlevs,&
 svc_tau,svc_lscale,sppt_tau,sppt_lscale,&
 iseed_svc,iseed_sppt,iseed_shum,shum,shum_tau,shum_lscale
 use patterngenerator, only: random_pattern, patterngenerator_init,&
 getnoise, patterngenerator_advance, patterngenerator_destroy
 use physcons, only: pi => con_pi
 use pressure_data, only: sl
 implicit none
 private
 public :: init_stochdata,destroy_stochdata
 complex(r_kind), allocatable, public, dimension(:) :: &
 spec_svc,spec_sppt,spec_shum
 real(r_kind), allocatable, public, dimension(:,:) :: &
 grd_svc,grd_sppt,grd_shum
 real(r_kind), allocatable, public, dimension(:) :: &
 vfact_shum,vfact_svc,vfact_sppt
 type(random_pattern), public :: rpattern_svc
 type(random_pattern), public :: rpattern_sppt
 type(random_pattern), public :: rpattern_shum

 contains
 subroutine init_stochdata()
    real(r_kind) delt,sigtop,sigbot
    integer n,nspinup,k
    delt = dt
    if (svc > tiny(svc)) then
       allocate(spec_svc(ndimspec))
       allocate(grd_svc(nlons,nlats))
       allocate(vfact_svc(nlevs))
       call patterngenerator_init(svc_lscale,delt,svc_tau,svc,iseed_svc,rpattern_svc,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = 10.*svc_tau/delt
       call getnoise(rpattern_svc,spec_svc)
       spec_svc = spec_svc*rpattern_svc%varspectrum
       do n=1,nspinup
          call patterngenerator_advance(spec_svc,rpattern_svc)
       enddo
       sigbot = 0.05; sigtop = 0.01
       do k=1,nlevs
          vfact_svc(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigbot) then
            vfact_svc(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_svc(k) = 0.
          endif
       enddo
    endif
    if (sppt > tiny(sppt)) then
       allocate(spec_sppt(ndimspec))
       allocate(grd_sppt(nlons,nlats))
       allocate(vfact_sppt(nlevs))
       call patterngenerator_init(sppt_lscale,delt,sppt_tau,sppt,iseed_sppt,rpattern_sppt,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = 10.*sppt_tau/delt
       call getnoise(rpattern_sppt,spec_sppt)
       spec_sppt = spec_sppt*rpattern_sppt%varspectrum
       do n=1,nspinup
          call patterngenerator_advance(spec_sppt,rpattern_sppt)
       enddo
       sigbot = 0.05; sigtop = 0.01
       do k=1,nlevs
          vfact_sppt(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigbot) then
            vfact_sppt(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_sppt(k) = 0.
          endif
       enddo
    endif
    if (shum > tiny(shum)) then
       allocate(spec_shum(ndimspec))
       allocate(grd_shum(nlons,nlats))
       allocate(vfact_shum(nlevs))
       call patterngenerator_init(shum_lscale,delt,shum_tau,shum,iseed_shum,rpattern_shum,&
                                  nlons,nlats,ntrunc,2.*pi)
       nspinup = 10.*shum_tau/delt
       call getnoise(rpattern_shum,spec_shum)
       spec_shum = spec_shum*rpattern_shum%varspectrum
       sigbot = 0.2
       ! humidity pert decays exponentially away from sfc
       do k=1,nlevs
          vfact_shum(k) = exp((sl(k)-1.)/sigbot)
       enddo
       do n=1,nspinup
          call patterngenerator_advance(spec_shum,rpattern_shum)
       enddo
    endif
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
 end subroutine destroy_stochdata

end module stoch_data
