      SUBROUTINE SFC_DIAG(IM,KM,PS,U1,V1,T1,Q1,
     &                  TSKIN,QSURF,F10M,U10M,V10M,T2M,Q2M,
     &                  PRSLKI,SLIMSK,EVAP,FM,FH,FM10,FH2)
!
      USE MACHINE , ONLY : kind_phys
      USE FUNCPHYS, ONLY : fpvs
      USE PHYSCONS, grav => con_g, SBC => con_sbc, HVAP => con_HVAP
     &,             CP => con_CP, HFUS => con_HFUS, JCAL => con_JCAL
     &,             EPS => con_eps, EPSM1 => con_epsm1
     &,             RVRDM1 => con_FVirt, RD => con_RD
      implicit none
!
      integer              IM, km
!
      real(kind=kind_phys) PS(IM),       U1(IM),      V1(IM),
     &                     T1(IM),       Q1(IM),  
     &                     TSKIN(IM),    QSURF(IM), 
     &                     F10M(IM),     U10M(IM),
     &                     V10M(IM),     T2M(IM),     Q2M(IM),
     &                                   PRSL1(IM),   PRSLKI(IM),
     &                     SLIMSK(IM),   EVAP(IM),    
     &                     FM(IM),       FH(IM),
     &                     FM10(IM),     FH2(IM)
!
!     Locals
!
      real (kind=kind_phys), parameter :: qmin=1.0e-8
      integer              k,i
!
      real(kind=kind_phys) QSS(IM), THETA1(IM)
!
      real(kind=kind_phys) g,    sig2k, fhi
!
      PARAMETER (G=grav)
!
      LOGICAL FLAG(IM), FLAGSNW(IM)
      real(kind=kind_phys) KT1(IM),       KT2(IM),      KTSOIL,
     &                     ET(IM,KM),
     &                     STSOIL(IM,KM), AI(IM,KM),    BI(IM,KM),
     &                     CI(IM,KM),     RHSTC(IM,KM)
!
!
!     ESTIMATE SIGMA ** K AT 2 M
!
      SIG2K = 1. - 4. * G * 2. / (CP * 280.)
!
!  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIE
!  PS IS IN PASCALS
!  THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1
!
!!
      DO I=1,IM
        THETA1(I) = T1(I) * PRSLKI(I)
      ENDDO
!!
!
      DO I = 1, IM
        F10M(I) = FM10(I) / FM(I)
        F10M(I) = min(F10M(I),1.)
        U10M(I) = F10M(I) *  U1(I)
        V10M(I) = F10M(I) * V1(I)
        FHI     = FH2(I) / FH(I)
        T2M(I)  = TSKIN(I)*(1. - FHI) + THETA1(I)*FHI
        T2M(I)  = T2M(I) * SIG2K
        IF(EVAP(I) >= 0.) THEN !  For EVAPORATION>0, USE INFERRED QSURF TO DEDUCE Q2M
          Q2M(I) = QSURF(I)*(1.-FHI) + max(qmin,Q1(I))*FHI
        ELSE                   !  FOR DEW FORMATION, USE SATURATED Q AT TSKIN
          qss(I) = fpvs(tskin(I))
          QSS(I) = EPS * QSS(I) / (PS(I) + EPSM1 * QSS(I))
          Q2M(I) = QSS(I)*(1.-FHI) + max(qmin,Q1(I))*FHI
        ENDIF
        QSS(I) = fpvs(t2m(I))
        QSS(I) = EPS * QSS(I) / (PS(I) + EPSM1 * QSS(I))
        Q2M(I) = MIN(Q2M(I),QSS(I))
      ENDDO

      RETURN
      END
