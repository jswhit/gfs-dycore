module params
! holds model parameters
! Public subroutines:
! read_namelist: read namelist
 use kinds, only: r_kind
 use sigio_module, only: sigio_head, sigio_srhead, sigio_sropen, sigio_sclose

 implicit none
 private

 public :: read_namelist,initfile,fhmax,dt,ntmax,ndimspec,nlons,nlats,&
 tstart,ndiss,efold,nlevs,ntrunc,sighead,dry,explicit,heldsuarez,jablowill,&
 ntout,fhout,idate_start

 character(len=500) :: initfile ! init cond filename
 integer            :: fhmax ! hours to run
 integer            :: fhout ! interval for IO
 real(r_kind)     :: dt    ! time step (secs)
 integer    :: ntmax ! time steps to run
 integer    :: nlons ! number of longitudes on grid
 integer    :: nlats ! number of latitudes on grid
 integer    :: nlevs ! number of levels on grid
 integer    :: ntrunc ! spectral truncation
 integer    :: ndimspec ! spectral array dimension
 type(sigio_head),save  :: sighead ! header struct from initfile
 logical    :: dry = .false.
 ! held-suarez forcing
 logical    :: heldsuarez = .false.
 ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
 ! idealized baroclinic instability test case.
 logical    :: jablowill = .false.
 ! use explicit time differencing
 ! if .true., explicit RK3 is used.
 ! if .false., semi-implicit RK3 (Kar, 2006, MWR p. 2916,
 ! http://dx.doi.org/10.1175/MWR3214.1) is used.
 ! explicit or semi-implicit rk3
 !logical    :: explicit = .true. ! use explicit rk3
 logical    :: explicit = .false. ! use semi-implicit rk3
 ! starting forecast time in seconds (read in from initfile)
 real(r_kind) :: tstart
 integer    :: idate_start(4) ! starting date (hr,month,day,year)
 integer    :: ntout ! time step interval for IO
 integer    :: ndiss=6 ! hyperdiffusion order
 real(r_kind) :: efold=3.*3600. ! efolding scale for smallest resolvable wave

 namelist/nam_dyn/initfile,fhmax,dt,dry,efold,ndiss,jablowill,heldsuarez,explicit,&
 fhout

 contains

 subroutine read_namelist()
   integer lu,iret
   initfile=""
   fhmax = 0
   fhout = 0
   dt = 0
   read(5,nam_dyn)
   if (initfile == "") then
      print *,'initfile must be specified'
      stop
   endif
   if (fhmax == 0) then
      print *,'fhmax must be specified'
      stop
   endif
   if (fhmax == 0) then
      print *,'fhout must be specified'
      stop
   endif
   if (dt == 0) then
      print *,'dt must be specified'
      stop
   endif
   if (mod(3600.d0,dt) .ne. 0) then
      print *,'1 hour must be an integer number of timesteps'
      stop
   endif
   ntmax = nint(real(fhmax)*3600./dt)
   ntout = nint(real(fhout)*3600/dt)
   lu = 7
   call sigio_sropen(lu,trim(initfile),iret)
   if (iret .ne. 0) then
      print *,'error opening ',trim(initfile),iret
      stop
   endif
   call sigio_srhead(lu,sighead,iret)
   if (iret .ne. 0) then
      print *,'error reading header from ',trim(initfile),iret
      stop
   else
      nlons = sighead%lonb
      nlats = sighead%latb
      nlevs = sighead%levs
      ntrunc = sighead%jcap
      tstart = sighead%fhour*3600.
      idate_start = sighead%idate
      print *,'nlons,nlats,nlevs,ntrunc=',nlons,nlats,nlevs,ntrunc
      print *,'tstart=',tstart,' secs'
      print *,'idate_start=',idate_start
   endif 
   if (jablowill .and. heldsuarez) then
      print *,'conflicting namelist options'
      print *,'heldsuarez and jablowill both cannot be .true.'
   endif
   ! for these idealized tests, model is dry.
   if (jablowill .or. heldsuarez) dry = .true.
   if (jablowill) print *,'running jablonowsky and williamson test..'
   if (heldsuarez) print *,'running held-suarez test..'
   call sigio_sclose(lu,iret)
   ndimspec = (ntrunc+1)*(ntrunc+2)/2
   print *,'namelist nam_dyn:'
   write(6, nam_dyn)
 end subroutine read_namelist

end module params
