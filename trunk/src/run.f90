module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).
use kinds, only: r_kind
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  ntout, explicit
use dyn_run, only: getdyntend, semimpadj
use phy_run, only: getphytend
use dyn_init, only: wrtout
use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
                          spfhumspec
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
  real(r_kind) t,t1,t2,fh
  real(r_kind), dimension(nlons,nlats,nlevs) :: spd
  integer count_0, count_1, count_rate, count_max
  character(len=500) filename

!$omp parallel
  my_id = omp_get_thread_num()
  if (my_id .eq. 0) print *,'running with',omp_get_num_threads(),' threads'
!$omp end parallel
  ! time step loop
  do nt=1,ntmax
     t = tstart + nt*dt
     call system_clock(count_0, count_rate, count_max)
     t1 = count_0*1.d0/count_rate
     ! advance solution with RK3
     call advance()
     call system_clock(count_1, count_rate, count_max)
     t2 = count_1*1.d0/count_rate
     spd = sqrt(ug**2+vg**2) ! max wind speed
     ! write out data at specified intervals.
     if (ntout .ne. 0 .and. mod(nt,ntout) .eq. 0) then
        fh = t/3600.
        write(filename,9000) int(fh)
9000    format('sig.f',i0.3) ! at least three digits used
        print *,'writing to ',trim(filename)
        call wrtout(t/3600.,filename)
     end if
     write(6,9001) t/3600.,maxval(spd),minval(psg/100.),maxval(psg/100.),t2-t1
9001 format('t = ',f8.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', cpu time = ',f7.3)
  enddo

  return
end subroutine run

subroutine advance()
! advance model state to next time step.
! (using explicit or semi-implicit third-order runge-kutta)
! hybrid sigma-pressure dynamical core described in
! http://www.emc.ncep.noaa.gov/officenotes/newernotes/on462.pdf
! The only difference is that here we use a forward in time
! semi-implicit runge-kutta scheme described by
! Kar (2006, http://journals.ametsoc.org/doi/pdf/10.1175/MWR3214.1)
! instead of semi-lmplicit assellin-filtered leap-frog.
  complex(r_kind),dimension(ndimspec,nlevs) :: &
  vrtspec_save,divspec_save,virtempspec_save,spfhumspec_save
  complex(r_kind),dimension(ndimspec) :: lnpsspec_save
  complex(r_kind), dimension(ndimspec,nlevs) :: &
  dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt
  complex(r_kind), dimension(ndimspec) :: dlnpsspecdt
  integer k
  real(r_kind) dtx
  ! save original fields.
  vrtspec_save = vrtspec
  divspec_save = divspec
  virtempspec_save = virtempspec
  spfhumspec_save = spfhumspec
  lnpsspec_save = lnpsspec
  do k=0,2
     dtx = dt/float(3-k)
     ! dynamics tendencies.
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt,dlnpsspecdt)
     ! add physics tendencies.
     call getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dspfhumspecdt,dlnpsspecdt)
     if (.not. explicit) then
         ! semi-implicit adjustment.
         call semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                        divspec_save,virtempspec_save,lnpsspec_save,k,dtx)
     endif
     vrtspec=vrtspec_save+dtx*dvrtspecdt
     divspec=divspec_save+dtx*ddivspecdt
     virtempspec=virtempspec_save+dtx*dvirtempspecdt
     spfhumspec=spfhumspec_save+dtx*dspfhumspecdt
     lnpsspec=lnpsspec_save+dtx*dlnpsspecdt
  enddo
end subroutine advance

end module run_mod
