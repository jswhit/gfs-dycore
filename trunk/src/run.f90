module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).
use kinds, only: r_kind
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  ntout, explicit, idate_start, adiabatic, ntrac
use dyn_run, only: getdyntend, semimpadj
use phy_run, only: getphytend
use dyn_init, only: wrtout
use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
                          tracerspec
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
  real(r_kind) t,fh
  real(8) t1,t2
  real(r_kind), dimension(nlons,nlats,nlevs) :: spd
  integer(8) count, count_rate, count_max
  character(len=500) filename

!$omp parallel
  my_id = omp_get_thread_num()
  if (my_id .eq. 0) print *,'running with',omp_get_num_threads(),' threads'
!$omp end parallel
  ! time step loop
  do nt=1,ntmax
     t = tstart + nt*dt
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     ! advance solution with RK3
     call advance()
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
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
9001 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', cpu time = ',f0.3)
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
  vrtspec_save,divspec_save,virtempspec_save
  complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: &
  tracerspec_save
  complex(r_kind),dimension(ndimspec) :: lnpsspec_save
  complex(r_kind), dimension(ndimspec,nlevs) :: &
  dvrtspecdt,ddivspecdt,dvirtempspecdt
  complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: &
  dtracerspecdt
  complex(r_kind), dimension(ndimspec) :: dlnpsspecdt
  integer k
  real(r_kind) dtx
  ! save original fields.
  vrtspec_save = vrtspec
  divspec_save = divspec
  virtempspec_save = virtempspec
  if (ntrac > 0) tracerspec_save = tracerspec
  lnpsspec_save = lnpsspec
  do k=0,2
     dtx = dt/float(3-k)
     ! dynamics tendencies.
     call getdyntend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt)
     ! add physics tendencies.
     if (.not. adiabatic) &
     call getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,dtx)
     if (.not. explicit) then
         ! semi-implicit adjustment.
         call semimpadj(ddivspecdt,dvirtempspecdt,dlnpsspecdt,&
                        divspec_save,virtempspec_save,lnpsspec_save,k,dtx)
     endif
     vrtspec=vrtspec_save+dtx*dvrtspecdt
     divspec=divspec_save+dtx*ddivspecdt
     virtempspec=virtempspec_save+dtx*dvirtempspecdt
     if (ntrac > 0) tracerspec=tracerspec_save+dtx*dtracerspecdt
     lnpsspec=lnpsspec_save+dtx*dlnpsspecdt
  enddo
end subroutine advance

!subroutine getvaliddate(tstart,fhour,idate_start,idate_valid)
!   ! Compute valid time from initial date and forecast hour
!   ! (using NCEP w3lib)
!   real(r_kind), intent(in) :: tstart
!   real(r_kind), intent(in) :: fhour
!   integer,dimension(8):: id,jd
!   real(r_kind), dimension(5):: fh
!   integer, intent(in),  dimension(4) :: idate_start
!   integer, intent(out), dimension(4) :: idate_valid
!   fh=zero; id=0; jd=0
!   fh(2)=fhour    ! relative time interval in hours
!   id(1)=idate_start(4) ! year
!   id(2)=idate_start(2) ! month
!   id(3)=idate_start(3) ! day
!   id(4)=0        ! time zone
!   id(5)=idate_start(1) ! hour
!   call w3movdat(fh,id,jd)
!   !     JDAT       INTEGER NCEP ABSOLUTE DATE AND TIME
!   !                (YEAR, MONTH, DAY, TIME ZONE,
!   !                 HOUR, MINUTE, SECOND, MILLISECOND)
!   idate_valid(1)=jd(5) ! hour
!   idate_valid(2)=jd(2) ! mon
!   idate_valid(3)=jd(3) ! day
!   idate_valid(4)=jd(1) ! year
!   return
!end subroutine getvaliddate

end module run_mod
