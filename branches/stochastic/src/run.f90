module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).
use kinds, only: r_kind,r_double
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  heldsuarez,jablowill,fhzer,ntrac,ntout, explicit, idate_start, adiabatic, ntrac,&
  sfcinitfile, postphys, ntdfi, svc, sppt, spdt, sppt_logit, spdt_logit
use shtns, only: lats, gauwts, spectogrd
use dyn_run, only: getdyntend, semimpadj
use phy_run, only: getphytend
use phy_data, only: wrtout_sfc, wrtout_flx, init_phydata
use dyn_init, only: wrtout_sig
use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
                          tracerspec, disspec, dmp_prof, diff_prof
! these arrays used to print diagnostics after each time step.
use pressure_data, only: psg
use grid_data, only: ug,vg,dlnpsdt
private
public :: run
contains

subroutine run()
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  implicit none
  integer i,nt,ntstart,my_id
  real(r_kind) fh,fha,pstendmean
  real(r_double) ta,t
  real(8) t1,t2
  real(r_kind), dimension(nlons,nlats,nlevs) :: spd
  real(r_kind), dimension(nlons,nlats) :: pstend,areawts
  complex(r_kind), dimension(:,:), allocatable :: &
  vrtspec_dfi,divspec_dfi,virtempspec_dfi
  complex(r_kind), dimension(:,:,:), allocatable :: tracerspec_dfi
  complex(r_kind), dimension(:), allocatable :: lnpsspec_dfi
  real(r_kind), dimension(:), allocatable :: dfi_wts
  integer(8) count, count_rate, count_max
  character(len=500) filename,filename_save

  do i=1,nlons
     areawts(i,:) = gauwts(:)
  enddo
  areawts = areawts/sum(areawts)

!$omp parallel
  my_id = omp_get_thread_num()
  if (my_id .eq. 0) print *,'running with',omp_get_num_threads(),' threads'
!$omp end parallel

  t = tstart
  ta = 0.
  ntstart = 1

  ! digital filter loop.
  if (ntdfi > 0) then
     print *,'in dfi time step loop...'
     ! allocate work space
     allocate(vrtspec_dfi(ndimspec,nlevs),divspec_dfi(ndimspec,nlevs))
     allocate(virtempspec_dfi(ndimspec,nlevs),lnpsspec_dfi(ndimspec))
     allocate(tracerspec_dfi(ndimspec,nlevs,ntrac),dfi_wts(0:2*ntdfi))
     ! compute dfi weights
     call set_dfi_wts(dfi_wts)
     ! intialize weighted time averages.
     vrtspec_dfi = dfi_wts(0)*vrtspec
     divspec_dfi = dfi_wts(0)*divspec
     virtempspec_dfi = dfi_wts(0)*virtempspec
     tracerspec_dfi = dfi_wts(0)*tracerspec
     lnpsspec_dfi = dfi_wts(0)*lnpsspec
     do nt=1,2*ntdfi
        call system_clock(count, count_rate, count_max)
        t1 = count*1.d0/count_rate
        call advance(t)
        t = t + dt ! absolute forecast time.
        ta = ta + dt ! absolute forecast time.
        vrtspec_dfi = vrtspec_dfi + dfi_wts(nt)*vrtspec
        divspec_dfi = divspec_dfi + dfi_wts(nt)*divspec
        virtempspec_dfi = virtempspec_dfi + dfi_wts(nt)*virtempspec
        tracerspec_dfi = tracerspec_dfi + dfi_wts(nt)*tracerspec
        lnpsspec_dfi = lnpsspec_dfi + dfi_wts(nt)*lnpsspec
        call system_clock(count, count_rate, count_max)
        t2 = count*1.d0/count_rate
        spd = sqrt(ug**2+vg**2) ! max wind speed
        pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
        pstendmean = sqrt(sum(pstend*areawts))
        write(6,9002) t/3600.,maxval(spd),minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
        ! write out surface and flux data in middle of dfi window.
        if (nt .eq. ntdfi) then
           fh = t/3600.
           write(filename_save,9000) nint(fh)
           print *,'writing to ',trim(filename_save),' fh=',fh
           call wrtout_sfc(fh,filename_save)
           write(filename,9001) nint(fh)
           print *,'writing to ',trim(filename),' fh=',fh
           call wrtout_flx(fh,ta,filename)
        ! write first time step output
        else if (nt .eq. 1) then
           write(filename,8999) nint(fh)
           print *,'writing to ',trim(filename),' fh=',fh
           call wrtout_sig(fh,filename)
           write(filename,9000) nint(fh)
           print *,'writing to ',trim(filename),' fh=',fh
           call wrtout_sfc(fh,filename)
           write(filename,9001) nint(fh)
           print *,'writing to ',trim(filename),' fh=',fh
           call wrtout_flx(fh,ta,filename)
        end if
     enddo
     print *,'done with dfi loop, resetting fields and restarting...'
     ! reset model state to weighted time average.
     vrtspec = vrtspec_dfi; divspec = divspec_dfi; virtempspec = virtempspec_dfi
     lnpsspec = lnpsspec_dfi; tracerspec = tracerspec_dfi
     ! deallocate work space.
     deallocate(vrtspec_dfi,divspec_dfi,virtempspec_dfi,lnpsspec_dfi,tracerspec_dfi,dfi_wts)
     ! reset surface data to values at middle of window (also zeros flux arrays).
     sfcinitfile = filename_save; call init_phydata(); ta = 0.
     ! reset time.
     t = tstart + ntdfi*dt; ntstart = ntdfi+1
     ! write out spectral data after dfi.
     fh = t/3600.
     write(filename,8999) nint(fh)
     print *,'writing to ',trim(filename),' fh=',fh
     call wrtout_sig(fh,filename)
  endif

  ! main time step loop
  do nt=ntstart,ntmax
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     ! advance solution with RK3
     call advance(t)
     t = t + dt ! absolute forecast time.
     ta = ta + dt ! time in accumulaton interval.
     fh = t/3600.
     fha = ta/3600.
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     spd = sqrt(ug**2+vg**2) ! max wind speed
     pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
     pstendmean = sqrt(sum(pstend*areawts))
     write(6,9002) t/3600.,maxval(spd),minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
9002 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', pstend = ',f0.3,', cpu time = ',f0.3)
     ! write out data at specified intervals.
     ! data always written at first time step.
     if (nt .eq. 1 .or. (ntout .ne. 0 .and. mod(nt,ntout) .eq. 0)) then
        write(filename,8999) nint(fh)
8999    format('SIG.F',i0.2) ! at least three digits used
        print *,'writing to ',trim(filename),' fh=',fh
        call wrtout_sig(fh,filename)
        ! write out boundary and flux files if using gfs physics.
        write(filename,9000) nint(fh)
9000    format('SFC.F',i0.2) ! at least three digits used
        print *,'writing to ',trim(filename),' fh=',fh
        call wrtout_sfc(fh,filename)
        write(filename,9001) nint(fh)
9001    format('FLX.F',i0.2) ! at least three digits used
        print *,'writing to ',trim(filename),' fh=',fh
        call wrtout_flx(fh,ta,filename)
     end if
     if (abs(fha-fhzer) .lt. 1.e-5) ta=0. ! reset accum time.
  enddo

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
  use patterngenerator, only: patterngenerator_advance
  use stoch_data, only: rpattern_svc,rpattern_sppt,rpattern_spdt,&
  spec_svc,spec_sppt,spec_spdt,grd_svc,grd_sppt,grd_spdt
  real(r_double), intent(in) :: t
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
  real(r_double) dtx
  logical :: profile = .true. ! print out timing stats
  integer(8) count, count_rate, count_max
  real(8) t1,t2
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
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     if (k .eq. 0) then
         if (svc > tiny(svc)) then
            call patterngenerator_advance(spec_svc,rpattern_svc)
            call spectogrd(spec_svc,grd_svc)
         endif
         if (sppt > tiny(sppt)) then
            call patterngenerator_advance(spec_sppt,rpattern_sppt)
            call spectogrd(spec_sppt,grd_sppt)
            ! logit transform to bounded interval [-1,+1]
            if (sppt_logit) grd_sppt = (2./(1.+exp(grd_sppt)))-1.
         endif
         if (spdt > tiny(spdt)) then
            call patterngenerator_advance(spec_spdt,rpattern_spdt)
            call spectogrd(spec_spdt,grd_spdt)
            ! logit transform to bounded interval [-1,+1]
            if (spdt_logit) grd_spdt = (2./(1.+exp(grd_spdt)))-1.
         endif
     endif
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt)
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     if (profile) print *,'time in getdyntend=',t2-t1
     ! compute physics tendencies at beginning, hold constant within RK3 sub time steps.
     ! (process-split approach)
     if (.not. adiabatic .and. .not. postphys) then
        if (k == 0) then
           call system_clock(count, count_rate, count_max)
           t1 = count*1.d0/count_rate
           call getphytend(dvrtspecdt_phy,ddivspecdt_phy,dvirtempspecdt_phy,dtracerspecdt_phy,dlnpsspecdt_phy,t,dt)
           call system_clock(count, count_rate, count_max)
           t2 = count*1.d0/count_rate
           if (profile) print *,'time in getphytend=',t2-t1
        endif
        dvrtspecdt = dvrtspecdt + dvrtspecdt_phy 
        ddivspecdt = ddivspecdt + ddivspecdt_phy
        dvirtempspecdt = dvirtempspecdt + dvirtempspecdt_phy
        dtracerspecdt = dtracerspecdt + dtracerspecdt_phy
        dlnpsspecdt = dlnpsspecdt + dlnpsspecdt_phy
     endif 
     ! semi-implicit adjustment includes physics tendencies for process-split
     ! physics.
     if (.not. explicit) then
         ! semi-implicit adjustment.
         call system_clock(count, count_rate, count_max)
         t1 = count*1.d0/count_rate
         call semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                        divspec_save,virtempspec_save,lnpsspec_save,k,dtx)
         call system_clock(count, count_rate, count_max)
         t2 = count*1.d0/count_rate
         if (profile) print *,'time in semimpadj=',t2-t1
     endif
     ! update for RK3 sub-step.
     vrtspec=vrtspec_save+dtx*dvrtspecdt
     divspec=divspec_save+dtx*ddivspecdt
     virtempspec=virtempspec_save+dtx*dvirtempspecdt
     if (ntrac > 0) tracerspec=tracerspec_save+dtx*dtracerspecdt
     lnpsspec=lnpsspec_save+dtx*dlnpsspecdt
     ! forward implicit treatment of linear damping/diffusion
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
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
      call system_clock(count, count_rate, count_max)
      t2 = count*1.d0/count_rate
      if (profile) print *,'time in diffusion update=',t2-t1
  enddo
  ! apply physics parameterizations as an adjustment to fields updated by dynamics.
  ! (time-split approach)
  if (.not. adiabatic .and. postphys) then
     ! update variables on grid.
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,&
          dtracerspecdt,dlnpsspecdt_phy,.true.) ! <- .true. for spectogrd only
     ! compute physics tendencies, apply as an adjustment to updated state
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     call getphytend(dvrtspecdt_phy,ddivspecdt_phy,dvirtempspecdt_phy,dtracerspecdt_phy,dlnpsspecdt_phy,t,dt)
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     if (profile) print *,'time in getphytend=',t2-t1
     vrtspec=vrtspec+dt*dvrtspecdt_phy
     divspec=divspec+dt*ddivspecdt_phy
     virtempspec=virtempspec+dt*dvirtempspecdt_phy
     if (ntrac > 0) tracerspec=tracerspec+dt*dtracerspecdt_phy
     !lnpsspec=lnpsspec+dt*dlnpsspecdt_phy ! physics does not modify ps
  end if
end subroutine advance

subroutine set_dfi_wts(dfi_wts)
 ! set Lanczos filter weights for digital filter 
 real(r_kind), intent(out) :: dfi_wts(0:2*ntdfi)
 real(r_kind) totsum,sx,wx
 integer kstep
 dfi_wts = 1.
 totsum = 0.
 do kstep=0,2*ntdfi
    sx     = acos(-1.)*(kstep-ntdfi)/ntdfi
    wx     = acos(-1.)*(kstep-ntdfi)/(ntdfi+1)
    if (kstep .NE. ntdfi) then
       dfi_wts(kstep) = sin(wx)/wx*sin(sx)/sx
    endif
    totsum = totsum + dfi_wts(kstep)
 enddo
 dfi_wts = dfi_wts/totsum
 print *,'lanczos dfi wts:'
 do kstep=0,2*ntdfi
    print *,kstep-ntdfi,dfi_wts(kstep)
 enddo
end subroutine set_dfi_wts

end module run_mod
