module phy_data
! public subroutines:
! init_phydata: allocate and populate arrays.
! destroy_phydata: deallocate arrays.
 use kinds, only: r_kind,r_double
 use params, only: nlons,nlats,nlevs,ndimspec,sfcinitfile,nmtvr,ntoz,ntclw,num_p3d,num_p2d,&
      ntrunc,sighead,fhswr,fhlwr,idate_start,fhzer,dt,gloopb_filter
 use shtns, only: lats
 use physcons, only : tgice => con_tice, con_pi
 implicit none
 private
 public :: init_phydata, destroy_phydata, wrtout_flx, wrtout_sfc
 real(r_kind), public, allocatable, dimension(:,:) :: apcp, precip, pwat

 contains

 subroutine init_phydata()
    allocate(precip(nlons,nlats))
    allocate(apcp(nlons,nlats))
    allocate(pwat(nlons,nlats))
    apcp = 0 ! precip accumulated from start of run
 end subroutine init_phydata

 subroutine destroy_phydata()
    deallocate(precip,apcp,pwat)
 end subroutine destroy_phydata

 subroutine wrtout_sfc(fhour,filename)
   implicit none
   real(r_kind), intent(in) :: fhour
   character(len=120) filename
 end subroutine wrtout_sfc

 subroutine wrtout_flx(fhour,ta,filename)
   implicit none
   real(r_kind),intent(in) :: fhour
   real(r_double),intent(in) :: ta
   character(len=120),intent(in) :: filename
   integer, parameter :: noflx=7 
   real(r_kind) secswr,seclwr,zhour,fha,dtsw,dtlw
   integer,PARAMETER :: NFLD=25
   integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
   real(4) RTIMER(NFLD),rtime,rtimsw,rtimlw
   real(4) colat1
   integer,PARAMETER :: IPRS=1,ITEMP=11,IZNLW=33,IMERW=34,ISPHUM=51,IPWAT=54,&
            IAPCP=61,IPCPR=59,ISNOWD=65,ICLDF=71,ICCLDF=72,&
            ISLMSK=81,IZORL=83,IALBDO=84,ISOILM=144,ICEMSK=91,&
            ISIK=92,                               &
            ILHFLX=121,ISHFLX=122,IZWS=124,IMWS=125,IGHFLX=155,&
            IUSWFC=160,IDSWFC=161,IULWFC=162,IDLWFC=163,&
            INSWFC=164,INLWFC=165,&
            IDSWVB=166,IDSWVD=167,IDSWNB=168,IDSWND=169,&
            ITMX=15,ITMN=16,IRNOF=90,IEP=145,&
            IQMX=204,IQMN=205,&
            ICLDWK=146,IZGW=147,IMGW=148,IHPBL=221,&
            IDSWF=204,IDLWF=205,IUSWF=211,IULWF=212,ICPCPR=214,&
            IUVBF=200,IUVBFC=201,ISUNTM=191,&
            ICSUSW=160,ICSDSW=161,ICSULW=162,ICSDLW=163,&
            ISFC=1,ITOA=8,IELEV=105,&
            ISGLEV=109,IDBLS=111,I2DBLS=112,ICOLMN=200,&
            IBLBL=209,IBLTL=210,IBLLYR=211,&
            ILCBL=212,ILCTL=213,ILCLYR=214,&
            IMCBL=222,IMCTL=223,IMCLYR=224,&
            IHCBL=232,IHCTL=233,IHCLYR=234,&
            ICVBL=242,ICVTL=243,ICVLYR=244,&
            ISLC=160,ISNOD=66,&
            ISLO=222,ISBS=198,ISNC=238,ICMM=179,&
            ISNOHF=229,ISMCWLT=219,ISMCREF=220,&
            ICNP=223,&
            IVEG=87,IVTP=225,ISTP=224,IUST=253,IHGT=7,&
            IRST=140,ICHH=208,ISRF=235,IEVBS=199,&
            IEVCW=200,ITRAN=210,ISTC=86,&
            INST=10,IWIN=2,IAVG=3,IACC=4,&
            ICEN=7,IFHOUR=1,IFDAY=2
   LOGICAL(1) LBM(nlons*nlats)
   CHARACTER G(200+nlons*nlats*(16+1)/8)
   INTEGER   IPUR(NFLD),ITLR(NFLD),ICEN2,IGEN
   DATA      IPUR/IULWF , IUSWF , IUSWF , IDSWF ,  ICLDF,   IPRS,&
                  IPRS, ITEMP ,  ICLDF,   IPRS,   IPRS, ITEMP ,&
                  ICLDF,   IPRS,   IPRS, ITEMP ,  IUVBF, IUVBFC ,&
                  IDSWF, ICSULW, ICSUSW, ICSDLW, ICSUSW, ICSDSW,&
                  ICSULW /
   DATA      ITLR/ITOA  , ITOA  , ISFC  , ISFC  , IHCLYR, IHCTL ,&
                  IHCBL , IHCTL , IMCLYR, IMCTL , IMCBL , IMCTL ,&
                  ILCLYR, ILCTL , ILCBL , ILCTL , ISFC  , ISFC ,&
                  ITOA  ,  ITOA ,  ITOA ,  ISFC , ISFC  , ISFC,&
                  ISFC /
   INTEGER   IENS(5),IDS(255),i,j,n,k,k4,l
   real(4) wrkga(nlons*nlats)

   call BAOPENWT(NOFLX,filename,ierr)
   if (ierr .ne. 0) print *,' iostat after baclose of flxf ile ',ierr

   fha = ta/3600.
   zhour = fhour - fha ! last time accum arrays zeroed
   !print *,'fha,fhour,zhour=',fha,fhour,zhour
   if (fhswr .lt. dt) then
      dtsw = dt
   else
      dtsw  = 3600.0 * fhswr
   endif
   if (fhlwr .lt. dt) then
      dtlw = dt
   else
      dtlw  = 3600.0 * fhlwr
   endif
   SECSWR=MAX(ta,DTSW)
   SECLWR=MAX(ta,DTLW)
   ICEN2 = sighead%icen2; IGEN = sighead%igen

   CALL IDSDEF(1,IDS)
   ids(IQMX)   = 5
   ids(IQMN)   = 5
   ids(IUVBF)  = 2
   ids(IUVBFC) = 2
   ids(icemsk) = 3      ! ICE CONCENTRATION ()
   ids(isik)   = 2      ! ICE THICKNESS (M)
   ids(IZORL)  = 4
   ids(IHGT)   = 3
   ids(IVEG)   = 2
   ids(IUST)   = 3
   ids(ICHH)   = 4
   ids(ICMM)   = 4
   ids(ISRF)   = 5
   ids(ITEMP)  = 3
   ids(ISPHUM) = 6
   ids(IZNLW)  = 2
   ids(IMERW)  = 2
   ids(ISNC)   = 3
   ids(ISTC)   = 4
   ids(ISOILM) = 4
   ids(ISNOD)  = 6
   ids(ISNOWD) = 5
   ids(ICNP)   = 5
   ids(IPCPR)  = 6
   ids(ICPCPR) = 6
   ids(IRNOF)  = 5                                            
   ids(ISMCWLT)  = 4
   ids(ISMCREF)  = 4
   ILPDS   = 28
   IF(ICEN2.EQ.2) ILPDS=45
   IENS(1) = 1
   IENS(2) = sighead%iens(1)
   IENS(3) = sighead%iens(2)
   IENS(4) = 1
   IENS(5) = 255
   IYR     = IDATE_START(4)
   IMO     = IDATE_START(2)
   IDA     = IDATE_START(3)
   IHR     = IDATE_START(1)
   IFHR    = NINT(ZHOUR)
   ITHR     = NINT(FHOUR)
   IF(FHOUR.GT.ZHOUR) THEN
     RTIME = 1./(3600.*(FHOUR-ZHOUR))
   ELSE
     RTIME = 0.
   ENDIF
   IF(SECSWR.GT.0.) THEN
     RTIMSW = 1./SECSWR
   ELSE
     RTIMSW = 1.
   ENDIF
   IF(SECLWR.GT.0.) THEN
     RTIMLW = 1./SECLWR
   ELSE
     RTIMLW = 1.
   ENDIF
   RTIMER    = RTIMSW
   RTIMER(1) = RTIMLW
   RTIMER(20)=RTIMLW       ! CSULF_TOA
   RTIMER(22)=RTIMLW       ! CSDLF_SFC
   RTIMER(25)=RTIMLW       ! CSULF_SFC

   ! FIRST COLATITUDE OF GRID IF IDRT=4 (RADIANS)
   colat1 = 0.5*con_pi-lats(1,1)

   LBM = .true.

   call twodtooned(precip,wrkga)
   print *,'min/max prate',minval(wrkga),maxval(wrkga),rtime
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IPCPR,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,IFHR,0,INST,0,0,ICEN2,IDS(IPCPR),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '1)Precipitation rate (kg/m**2/s) land and sea surface        '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(apcp,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IAPCP,ISFC,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,0,ITHR,IACC,0,0,ICEN2,IDS(IAPCP),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '2)Total precip (kg/m**2) land sea surface '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   call twodtooned(pwat,wrkga)
   call gribit(wrkga,LBM,4,nlons,nlats,16,colat1,ILPDS,2,ICEN,IGEN,&
               0,IPWAT,ICOLMN,0,0,IYR,IMO,IDA,IHR,&
               IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPWAT),IENS,&
               0.,0.,0.,0.,0.,0.,G,LG,IERR)
   if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',&
    '3)Precipitable water (kg/m**2) total atmospheric column      '
   IF(IERR.EQ.0) CALL WRYTE(noflx,LG,G)

   PRINT *,'GRIB FLUX FILE WRITTEN ',fhour,idate_start
   call baclose(noflx,ierr)
   print *,' iostat after baclose of flxf file ',ierr

 end subroutine wrtout_flx

 subroutine twodtooned(data2,data1)
   real(r_kind), intent(in) :: data2(nlons,nlats)
   real(4), intent(out) :: data1(nlons*nlats)
   integer i,j,n
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      data1(n) = data2(i,j)
   enddo
 end subroutine twodtooned

end module phy_data
