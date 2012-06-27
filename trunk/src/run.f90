module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).
use kinds, only: r_kind
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  heldsuarez,jablowill,fhzer,ntrac,ntout, explicit, idate_start, adiabatic, ntrac
use dyn_run, only: getdyntend, semimpadj
use phy_run, only: getphytend
use phy_data, only: wrtout_sfc, wrtout_flx
use dyn_init, only: wrtout_sig
use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
                          tracerspec, disspec, dmp_prof, diff_prof
! these arrays used to print diagnostics after each time step.
use pressure_data, only: psg
use grid_data, only: ug,vg
private
public :: run
contains

subroutine run()
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  implicit none
  integer nt,my_id
  real(r_kind) t,fh,ta,fha
  real(8) t1,t2
  real(r_kind), dimension(nlons,nlats,nlevs) :: spd
  integer(8) count, count_rate, count_max
  character(len=500) filename

!$omp parallel
  my_id = omp_get_thread_num()
  if (my_id .eq. 0) print *,'running with',omp_get_num_threads(),' threads'
!$omp end parallel
  ! time step loop
  t = tstart
  ta = 0.
  do nt=1,ntmax
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     ! advance solution with RK3
     call advance(t)
     t = t + dt ! absolute forecast time.
     ta = ta + dt ! time in accumulaton interval.
     fh = t/3600.
     fha = ta/3600.
     if (abs(fh-fhzer) .lt. tiny(dt)) ta=0. ! reset accum time.
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     spd = sqrt(ug**2+vg**2) ! max wind speed
     ! write out data at specified intervals.
     ! data always written at first time step.
     if (nt .eq. 1 .or. (ntout .ne. 0 .and. mod(nt,ntout) .eq. 0)) then
        write(filename,8999) int(fh)
8999    format('sig.f',i0.3) ! at least three digits used
        print *,'writing to ',trim(filename)
        call wrtout_sig(t/3600.,filename)
        ! write out boundary and flux files if using gfs physics.
        if (.not. heldsuarez .and. .not. jablowill) then
        write(filename,9000) int(fh)
9000    format('sfc.f',i0.3) ! at least three digits used
        print *,'writing to ',trim(filename)
        call wrtout_sfc(t/3600.,filename)
        write(filename,9001) int(fh)
9001    format('flx.f',i0.3) ! at least three digits used
        print *,'writing to ',trim(filename)
        call wrtout_flx(t/3600.,ta,filename)
        endif
     end if
     write(6,9002) t/3600.,maxval(spd),minval(psg/100.),maxval(psg/100.),t2-t1
9002 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', cpu time = ',f0.3)
  enddo

  return
end subroutine run

subroutine advance(t)
! advance model state to next time step.
! (using explicit or semi-implicit third-order runge-kutta)
! hybrid sigma-pressure dynamical core described in
! http://www.emc.ncep.noaa.gov/officenotes/newernotes/on462.pdf
! The only difference is that here we use a forward in time
! semi-implicit runge-kutta scheme described by
! Kar (2006, http://journals.ametsoc.org/doi/pdf/10.1175/MWR3214.1)
! instead of semi-lmplicit assellin-filtered leap-frog.
  real(r_kind), intent(in) :: t
  complex(r_kind),dimension(ndimspec,nlevs) :: &
  vrtspec_save,divspec_save,virtempspec_save
  complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: &
  tracerspec_save
  complex(r_kind),dimension(ndimspec) :: lnpsspec_save
  complex(r_kind), dimension(ndimspec,nlevs) :: &
  dvrtspecdt,ddivspecdt,dvirtempspecdt,dvrtspecdt_phy,ddivspecdt_phy,dvirtempspecdt_phy
  complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: &
  dtracerspecdt,dtracerspecdt_phy
  complex(r_kind), dimension(ndimspec) :: dlnpsspecdt,dlnpsspecdt_phy
  integer k, nt, kk
  real(r_kind) dtx
  logical :: postphys = .false.
  ! save original fields.
  vrtspec_save = vrtspec
  divspec_save = divspec
  virtempspec_save = virtempspec
  if (ntrac > 0) tracerspec_save = tracerspec
  lnpsspec_save = lnpsspec
  ! update dynamics using RK3.
  do k=0,2
     dtx = dt/float(3-k)
     ! dynamics tendencies.
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt)
     if (.not. explicit) then
         ! semi-implicit adjustment.
         call semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                        divspec_save,virtempspec_save,lnpsspec_save,k,dtx)
     endif
     ! compute physics tendencies at beginning, hold constant within RK3 sub time steps.
     if (.not. adiabatic .and. .not. postphys) then
        if (k == 0) then
           call getphytend(dvrtspecdt_phy,ddivspecdt_phy,dvirtempspecdt_phy,dtracerspecdt_phy,dlnpsspecdt_phy,t,dt)
        endif
        dvrtspecdt = dvrtspecdt + dvrtspecdt_phy 
        ddivspecdt = ddivspecdt + ddivspecdt_phy
        dvirtempspecdt = dvirtempspecdt + dvirtempspecdt_phy
        dtracerspecdt = dtracerspecdt + dtracerspecdt_phy
        dlnpsspecdt = dlnpsspecdt + dlnpsspecdt_phy
     endif 
     ! update
     vrtspec=vrtspec_save+dtx*dvrtspecdt
     divspec=divspec_save+dtx*ddivspecdt
     virtempspec=virtempspec_save+dtx*dvirtempspecdt
     if (ntrac > 0) tracerspec=tracerspec_save+dtx*dtracerspecdt
     lnpsspec=lnpsspec_save+dtx*dlnpsspecdt
     ! forward implicit treatment of linear damping/diffusion
     !$omp parallel do private(kk,nt)
     do kk=1,nlevs
        vrtspec(:,kk) = vrtspec(:,kk)/(1. - (disspec(:)*diff_prof(kk) - dmp_prof(kk))*dtx)
        divspec(:,kk) = divspec(:,kk)/(1. - (disspec(:)*diff_prof(kk) - dmp_prof(kk))*dtx)
        virtempspec(:,kk) = virtempspec(:,kk)/(1. - disspec(:)*diff_prof(kk)*dtx)
        do nt=1,ntrac
           tracerspec(:,kk,nt) = tracerspec(:,kk,nt)/(1. - disspec(:)*diff_prof(kk)*dtx)
        enddo
      enddo
      !$omp end parallel do
  enddo
  ! apply physics parameterizations as an adjustment to fields updated by dynamics.
  if (.not. adiabatic .and. postphys) then
     ! update variables on grid.
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,&
          dtracerspecdt,dlnpsspecdt_phy,.true.) ! <- .true. for spectogrd only
     ! compute physics tendencies, apply as an adjustment to updated state
     call getphytend(dvrtspecdt_phy,ddivspecdt_phy,dvirtempspecdt_phy,dtracerspecdt_phy,dlnpsspecdt_phy,t,dt)
     vrtspec=vrtspec+dt*dvrtspecdt_phy
     divspec=divspec+dt*ddivspecdt_phy
     virtempspec=virtempspec+dt*dvirtempspecdt_phy
     if (ntrac > 0) tracerspec=tracerspec+dt*dtracerspecdt_phy
     !lnpsspec=lnpsspec+dt*dlnpsspecdt_phy ! physics does not modify ps
  end if
end subroutine advance

end module run_mod
