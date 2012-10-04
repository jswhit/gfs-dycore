 module phy_run
! compute physics tendencies for gfs physics
! getphytend: compute total physics tendencies.
! Public subroutines:
! getphytend: compute tendencies in spectral space.

 use params, only: nlevs,nlons,nlats,ntrunc,ndimspec,ntrac,dt,dry,dcmip,massfix,pdryini
 use kinds, only: r_kind,r_single,r_double
 use shtns, only: grdtospec, getvrtdivspec, lons, lats, areawts
 use grid_data, only: virtempg,dlnpdtg,tracerg,ug,vg
 use pressure_data, only:  prs,psg,pk,ak,bk,dpk
 use spectral_data, only: lnpsspec
 use phy_data, only: precip,apcp,pwat
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
               eps => con_eps, omega => con_omega, cvap => con_cvap, &
               grav => con_g, pi => con_pi, fv => con_fvirt, rk => con_rocp

 implicit none
 private

 public :: getphytend

 contains

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,t,dtx)
   use omp_lib, only: omp_get_num_threads, omp_get_thread_num
   ! compute physics tendencies for gfs physics
   integer, parameter :: r8 = selected_real_kind(12)
   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
   dvrtspecdt,ddivspecdt,dvirtempspecdt
   complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
   dtracerspecdt
   real(r_double), intent(in) :: t,dtx
   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
   real(r_kind), parameter :: qmin = 1.e-10 ! min value for clipping tracers
   real(r8) gt(nlevs),prsl(nlevs),prsi(nlevs+1),delp(nlevs),rdelp(nlevs),gq(nlevs),&
            dtp,lat(1),ps(1),precl(1),gu(nlevs),gv(nlevs)
   real(r_kind), allocatable, dimension(:,:,:) :: dtdt,dudt,dvdt
   real(r_kind), allocatable, dimension(:,:) :: wrkg
   complex(r_kind), allocatable, dimension(:) :: wrkspec
   real(4), allocatable, dimension(:,:,:) :: work4
   real(r_kind), allocatable, dimension(:,:,:,:) :: dtracersdt
   integer :: i,j,k,n,nt,testcase
   real(8) tstart,tend
   real(r_kind) :: pdry, pwatg, pcorr
   integer(8) count, count_rate, count_max
   logical :: testomp=.false.  ! openmp debug flag
   logical :: turbflux

   dtp = dtx

   if (dcmip/10 .eq. 4) testcase=1
   if (dcmip/10 .eq. 5) testcase=0

   turbflux = .false.  ! test case 4-1 or 4-2
   if (dcmip .eq. 43 .or. dcmip/10 .eq. 5) turbflux = .true.  ! cases 4-3 or 5.

   if (.not. dry) then

   ! simple moist physics and boundary layer mixing.

   allocate(dtdt(nlons,nlats,nlevs),dudt(nlons,nlats,nlevs),dvdt(nlons,nlats,nlevs))
   allocate(dtracersdt(nlons,nlats,nlevs,ntrac))

   call system_clock(count, count_rate, count_max)
   tstart = count*1.d0/count_rate
! physics loop over horiz. grid points.
!$omp parallel do private(n,k,i,j,&
!$omp& lat,ps,precl,gt,gu,gv,gq,prsi,prsl,delp,rdelp) schedule(dynamic)
   do n=1,nlons*nlats
      ! n=i+(j-1)*nlons
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      ! compute sensible temp (clip humidity in computation).
      do k=1,nlevs
         gq(nlevs-k+1) = tracerg(i,j,k,1)
         gu(nlevs-k+1) = ug(i,j,k)
         gv(nlevs-k+1) = vg(i,j,k)
         gt(nlevs-k+1) = virtempg(i,j,k)/(1.+fv*max(qmin,tracerg(i,j,k,1)))
         prsl(nlevs-k+1) = prs(i,j,k)
      enddo
      prsi = pk(i,j,:)
      do k=1,nlevs
         delp(k)=prsi(k+1)-prsi(k)
         rdelp(k) = 1./delp(k)
      enddo
      lat(1) = lats(i,j); ps(1) = psg(i,j)
      call simple_physics(1, nlevs, dtp, lat, gt, gq, gu, gv, prsl, prsi, &
                          delp, rdelp, ps, precl, testcase, turbflux)
      psg(i,j) = ps(1)
      precip(i,j) = precl(1)*1000. ! convert from m/s to mm/s
      ! convert sensible temp back to virt temp.
      ! (clip humidity in conversion)
      do k=1,nlevs
         gt(k) = gt(k)*(1.+fv*max(qmin,gq(k)))
      enddo
      ! compute tendencies, 
      ! update grid data.
      do k=1,nlevs
         dtdt(i,j,k) = gt(nlevs-k+1)-virtempg(i,j,k)
         dtracersdt(i,j,k,1) = gq(nlevs-k+1)-tracerg(i,j,k,1)
         dudt(i,j,k) = gu(nlevs-k+1)-ug(i,j,k)
         dvdt(i,j,k) = gv(nlevs-k+1)-vg(i,j,k)
      enddo
      do k=1,nlevs
         gt(k) = dpk(i,j,nlevs-k+1)
      enddo 
      pwat(i,j) = sum(gt*(dtracersdt(i,j,:,1)+tracerg(i,j,:,1)))
   enddo ! end loop over horiz grid points
!$omp end parallel do 
   pwat = pwat/grav
   print *,'min/max dtdt',minval(dtdt),maxval(dtdt)
   print *,'min/max dudt',minval(dudt),maxval(dudt)
   print *,'min/max dvdt',minval(dvdt),maxval(dvdt)
   print *,'min/max dtracer1dt',minval(dtracersdt(:,:,:,1)),maxval(dtracersdt(:,:,:,1))
   print *,'min/max tracer1',minval(tracerg(:,:,:,1)),maxval(tracerg(:,:,:,1))
   print *,'min/max pwat',minval(pwat),maxval(pwat)
   call system_clock(count, count_rate, count_max)
   tend = count*1.d0/count_rate
   print *,'time in simple physics = ',tend-tstart
   if (testomp) then
!$omp parallel
      i = omp_get_num_threads()
!$omp end parallel
      allocate(work4(nlons,nlats,nlevs))
      open(7,file='gbphys.dat',form='unformatted')
      if (i > 1) then
         read(7) work4
         print *,'dtdt diff',maxval(abs(work4-dtdt))
         read(7) work4
         print *,'dudt diff',maxval(abs(work4-dudt))
         read(7) work4
         print *,'dvdt diff',maxval(abs(work4-dvdt))
         do nt=1,ntrac
            read(7) work4
            !print *,minval(work4),maxval(work4)
            print *,nt,'dtracerdt diff',maxval(abs(work4-dtracersdt(:,:,:,nt)))
         enddo
      else
         work4 = dtdt
         write(7) work4
         work4 = dudt
         write(7) work4
         work4 = dvdt
         write(7) work4
         do nt=1,ntrac
            work4 = dtracersdt(:,:,:,nt)
            write(7) work4
         enddo
      endif
      close(7)
      deallocate(work4)
      stop
   endif
! compute physics tendencies in spectral space
   dtracerspecdt = 0.
!$omp parallel do private(k,nt)
   do k=1,nlevs
      call grdtospec(dtdt(:,:,k), dvirtempspecdt(:,k))
      call getvrtdivspec(dudt(:,:,k),dvdt(:,:,k),dvrtspecdt(:,k),ddivspecdt(:,k),rerth)
      call grdtospec(dtracersdt(:,:,k,1), dtracerspecdt(:,k,1))
   enddo
!$omp end parallel do 
   dlnpsspecdt=0 ! physics does not change surface pressure
   dvirtempspecdt = dvirtempspecdt/dtx
   dvrtspecdt = dvrtspecdt/dtx
   ddivspecdt = ddivspecdt/dtx
   dtracerspecdt = dtracerspecdt/dtx
   deallocate(dtdt,dudt,dvdt)
   deallocate(dtracersdt)
   apcp = apcp + precip*dtx
   pwatg = sum(areawts*pwat)
   print *,'min/max inst. precip rate (mm/sec) =',minval(precip),maxval(precip)
   print *,'min/max accum precip (mm) =',minval(apcp),maxval(apcp)
   print *,'global mean pwat (mm) =',pwatg

   else

   ! dry case
   dvirtempspecdt = 0.
   dvrtspecdt = 0.
   ddivspecdt = 0.
   dtracerspecdt = 0.
   dlnpsspecdt = 0.
   apcp = 0.
   precip = 0.
   pwat = 0.
   pwatg = 0.

   endif

! global mean dry mass 'fixer'
   if (massfix) then
! compute global mean dry ps.
      pdry = sum(areawts*psg) - grav*pwatg
      print *,'pdry after physics update',pdry
! implied ps correction needed to return dry mass to initial value
      pcorr = pdry - pdryini
! add constant correction to every grid point
      allocate(wrkg(nlons,nlats),wrkspec(ndimspec))
      wrkg = psg - pcorr 
! compute implied lnps tendency in spectral space.
      wrkg = log(wrkg)
      call grdtospec(wrkg,wrkspec)
      dlnpsspecdt = (wrkspec - lnpsspec)/dtx
      deallocate(wrkg,wrkspec)
   endif ! massfix

   return
 end subroutine getphytend

end module phy_run
