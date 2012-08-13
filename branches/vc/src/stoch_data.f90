module stoch_data

 use kinds, only: r_kind
 use params, only: dt,svc,sppt,spdt,ndimspec,nlons,nlats,nlevs,&
 svc_tau,svc_lscale,sppt_tau,sppt_lscale,spdt_tau,spdt_lscale,&
 iseed_svc,iseed_sppt,iseed_spdt
 use patterngenerator, only: random_pattern, patterngenerator_init,&
 getnoise, patterngenerator_advance
 use pressure_data, only: sl
 implicit none
 private
 public :: init_stochdata,destroy_stochdata
 complex(r_kind), allocatable, public, dimension(:) :: &
 spec_svc,spec_sppt,spec_spdt
 real(r_kind), allocatable, public, dimension(:,:) :: &
 grd_svc,grd_sppt,grd_spdt
 real(r_kind), allocatable, public, dimension(:) :: &
 vfact_svc,vfact_sppt,vfact_spdt
 type(random_pattern), public :: rpattern_svc
 type(random_pattern), public :: rpattern_sppt
 type(random_pattern), public :: rpattern_spdt

 contains
 subroutine init_stochdata()
    real(r_kind) delt,sigtop,sigbot
    integer n,nspinup,k
    dt = delt
    if (svc > 0.) then
       allocate(spec_svc(ndimspec))
       allocate(grd_svc(nlons,nlats))
       allocate(vfact_svc(nlevs))
       call patterngenerator_init(svc_lscale,delt,svc_tau,svc,iseed_svc,rpattern_svc)
       nspinup = 10.*svc_tau/delt
       call getnoise(spec_svc)
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
    if (sppt > 0.) then
       allocate(spec_sppt(ndimspec))
       allocate(grd_sppt(nlons,nlats))
       allocate(vfact_sppt(nlevs))
       call patterngenerator_init(sppt_lscale,delt,sppt_tau,sppt,iseed_sppt,rpattern_sppt)
       nspinup = 10.*sppt_tau/delt
       call getnoise(spec_sppt)
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
    if (spdt > 0.) then
       allocate(spec_spdt(ndimspec))
       allocate(grd_spdt(nlons,nlats))
       allocate(vfact_spdt(nlevs))
       call patterngenerator_init(spdt_lscale,delt,spdt_tau,spdt,iseed_spdt,rpattern_spdt)
       nspinup = 10.*spdt_tau/delt
       call getnoise(spec_spdt)
       spec_spdt = spec_spdt*rpattern_spdt%varspectrum
       do n=1,nspinup
          call patterngenerator_advance(spec_spdt,rpattern_spdt)
       enddo
       sigbot = 0.05; sigtop = 0.01
       do k=1,nlevs
          vfact_spdt(k) = 1.0
          if (sl(k) .lt. sigbot .and. sl(k) .gt. sigbot) then
            vfact_spdt(k) = (sl(k)-sigtop)/(sigbot-sigtop)
          else if (sl(k) .le. sigtop) then
            vfact_spdt(k) = 0.
          endif
       enddo
    endif
 end subroutine init_stochdata
 subroutine destroy_stochdata()
    if (allocated(spec_svc)) deallocate(spec_svc,vfact_svc)
    if (allocated(grd_svc)) deallocate(grd_svc)
    if (allocated(spec_sppt)) deallocate(spec_sppt,vfact_sppt)
    if (allocated(grd_sppt)) deallocate(grd_sppt)
    if (allocated(spec_spdt)) deallocate(spec_spdt,vfact_spdt)
    if (allocated(grd_spdt)) deallocate(grd_spdt)
 end subroutine destroy_stochdata

end module stoch_data
