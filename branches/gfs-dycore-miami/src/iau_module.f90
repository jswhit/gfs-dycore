module iau_module

 use kinds, only: r_kind, r_double
 use params, only: iaufiles_fg,iaufiles_anl,iaufhrs,iau,ndimspec,nlevs,ntrac,&
                   iau_delthrs
 use dyn_init, only: readin_sig

 private

 public :: init_iau, destroy_iau, getiauforcing

 complex(r_kind), dimension(:,:,:),allocatable ::  vrtspec,divspec,virtempspec
 complex(r_kind), dimension(:,:,:,:),allocatable :: tracerspec
 complex(r_kind), dimension(:,:), allocatable :: lnpsspec
 integer, public :: nfiles
 logical, public :: init_iauialized = .false.

 contains

 subroutine init_iau()
   use sigio_module, only: sigio_head
   complex(r_kind), dimension(ndimspec,nlevs) :: &
     vrtspec_tmp,divspec_tmp,virtempspec_tmp
   complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: tracerspec_tmp
   complex(r_kind), dimension(ndimspec) :: lnpsspec_tmp, topospec
   type(sigio_head) sighead
   integer, allocatable, dimension(:) :: idt
   character(len=120) filename
   integer n,nfilesall
   init_iauialized = .true.
   nfilesall = size(iaufiles_anl)
   nfiles = 0
   do n=1,nfilesall
      filename = iaufiles_anl(n)
      if (trim(filename) .eq. '' .or. iaufhrs(n) .lt. 0) exit
      nfiles = nfiles + 1
   enddo
   print *,'nfiles = ',nfiles
   if (nfiles < 2) then
     print *,'must be at least two files in iaufiles_fg and iaufiles_anal'
     stop
   endif
   allocate(idt(nfiles-1))
   idt = iaufhrs(2:nfiles)-iaufhrs(1:nfiles-1)
   do n=1,nfiles-1
      if (idt(n) .ne. iaufhrs(2)-iaufhrs(1)) then
        print *,'forecast intervals in iaufhrs must be constant'
        stop
      endif
   enddo
   print *,'iau interval = ',iau_delthrs,' hours'
   deallocate(idt)
   allocate(vrtspec(ndimspec,nlevs,nfiles))
   allocate(divspec(ndimspec,nlevs,nfiles))
   allocate(virtempspec(ndimspec,nlevs,nfiles))
   allocate(tracerspec(ndimspec,nlevs,ntrac,nfiles))
   allocate(lnpsspec(ndimspec,nfiles))
   do n=1,nfiles
      filename = iaufiles_fg(n)
      print *,'reading ',trim(filename)
      call readin_sig(trim(filename),vrtspec_tmp,divspec_tmp,virtempspec_tmp,&
                      tracerspec_tmp,lnpsspec_tmp,topospec,sighead)
      filename = iaufiles_anl(n)
      print *,'reading ',trim(filename)
      call readin_sig(trim(filename),vrtspec(:,:,n),divspec(:,:,n),virtempspec(:,:,n),&
                      tracerspec(:,:,:,n),lnpsspec(:,n),topospec,sighead)
!$omp workshare
      vrtspec(:,:,n) = vrtspec(:,:,n) - vrtspec_tmp
      divspec(:,:,n) = divspec(:,:,n) - divspec_tmp
      virtempspec(:,:,n) = virtempspec(:,:,n) - virtempspec_tmp
      tracerspec(:,:,:,n) = tracerspec(:,:,:,n) - tracerspec_tmp
      lnpsspec(:,n) = lnpsspec(:,n) - lnpsspec_tmp
!$omp end workshare
   enddo
 end subroutine init_iau

 subroutine getiauforcing(dvrtspecdt_iau,ddivspecdt_iau,dvirtempspecdt_iau,dtracerspecdt_iau,dlnpsspecdt_iau,t)
   complex(r_kind), dimension(ndimspec,nlevs), intent(out) :: &
     dvrtspecdt_iau,ddivspecdt_iau,dvirtempspecdt_iau
   complex(r_kind), dimension(ndimspec,nlevs,ntrac), intent(out) :: dtracerspecdt_iau
   complex(r_kind), dimension(ndimspec), intent(out) :: dlnpsspecdt_iau
   real(r_double), intent(in) :: t
   real(r_double) delt, dt
   integer n
   dvrtspecdt_iau = 0.; ddivspecdt_iau = 0.
   dvirtempspecdt_iau = 0; dlnpsspecdt_iau = 0.
   dtracerspecdt_iau = 0.
   dt = iau_delthrs*3600.
   if (t < iaufhrs(1)*3600. .or. t > iaufhrs(nfiles)*3600.) then
      print *,'no iau forcing'
      return
   endif
   if (t .eq. 3600.*iaufhrs(nfiles)) then
!$omp workshare
     dvrtspecdt_iau = vrtspec(:,:,nfiles)/dt
     ddivspecdt_iau = divspec(:,:,nfiles)/dt
     dvirtempspecdt_iau = virtempspec(:,:,nfiles)/dt
     dtracerspecdt_iau = tracerspec(:,:,:,nfiles)/dt
!$omp end workshare
     dlnpsspecdt_iau = lnpsspec(:,nfiles)/dt
     return
   else if (t .eq. 3600.*iaufhrs(1)) then
!$omp workshare
     dvrtspecdt_iau = vrtspec(:,:,1)/dt
     ddivspecdt_iau = divspec(:,:,1)/dt
     dvirtempspecdt_iau = virtempspec(:,:,1)/dt
     dtracerspecdt_iau = tracerspec(:,:,:,1)/dt
!$omp end workshare
     dlnpsspecdt_iau = lnpsspec(:,1)/dt
     return
   endif
   do n=1,nfiles
      if (iaufhrs(n)*3600. > t) exit
   enddo
   print *,'n,t,to',n,t/3600.,iaufhrs(n)
   delt = (iaufhrs(n)-(t/3600.))/(iaufhrs(n)-iaufhrs(n-1))
!$omp workshare
   dvrtspecdt_iau = ((1.-delt)*vrtspec(:,:,n) + delt*vrtspec(:,:,n-1))/dt
   ddivspecdt_iau = ((1.-delt)*divspec(:,:,n) + delt*divspec(:,:,n-1))/dt
   dvirtempspecdt_iau = ((1.-delt)*virtempspec(:,:,n) + delt*virtempspec(:,:,n-1))/dt
   dtracerspecdt_iau = ((1.-delt)*tracerspec(:,:,:,n) + delt*tracerspec(:,:,:,n-1))/dt
!$omp end workshare
   dlnpsspecdt_iau = ((1.-delt)*lnpsspec(:,n) + delt*lnpsspec(:,n-1))/dt
   print *,'getiauforcing:',t/3600.,1.-delt,n,iaufhrs(n),delt,n-1,iaufhrs(n-1)
 end subroutine getiauforcing

 subroutine destroy_iau()
   deallocate(vrtspec)
   deallocate(divspec)
   deallocate(virtempspec)
   deallocate(tracerspec)
   deallocate(lnpsspec)
 end subroutine destroy_iau

end module iau_module
