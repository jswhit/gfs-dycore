      subroutine vcnhyb(im,km,nm,dt,zint,zdot,zadv)
!                .      .    .                                       .
! subprogram:    vcnhyb      vertical advection instability filter
!   prgmmr: iredell          org: w/nmc23    date: 91-05-07
!
! abstract: filters vertical advection tendencies
!   in the dynamics tendency equation in order to ensure stability
!   when the vertical velocity exceeds the cfl criterion.
!   the vertical velocity in this case is sigmadot.
!   for simple second-order centered eulerian advection,
!   filtering is needed when vcn=zdot*dt/dz>1.
!   the maximum eigenvalue of the linear advection equation
!   with second-order implicit filtering on the tendencies
!   is less than one for all resolvable wavenumbers (i.e. stable)
!   if the nondimensional filter parameter is nu=(vcn**2-1)/4.
!
! program history log:
!   97-07-30  iredell
!
! usage:    call vcnhyb(im,km,nm,dt,zint,zdot,zadv,nvcn,xvcn)
!
!   input argument list:
!     im       - integer number of gridpoints to filter
!     km       - integer number of vertical levels
!     nm       - integer number of fields
!     dt       - real(r_kind) timestep in seconds
!     zint     - real(r_kind) (im,km+1) interface vertical coordinate values (top to bot)
!     zdot     - real(r_kind) (im,km+1) vertical coordinate velocity (top to bot)
!     zadv     - real(r_kind) (im,km,nm) vertical advection tendencies (bot to top)
!
!   output argument list:
!     zadv     - real(r_kind) (im,km,nm) vertical advection tendencies
!
! local vars:
!     nvcn     - integer number of points requiring filtering
!     xvcn     - real(r_kind) maximum vertical courant number
!
!   subprograms called:
!     tridim_hyb   - tridiagonal matrix solver
!
      use kinds, only: r_kind,r_double
      implicit none
      integer,intent(in):: im,km,nm
      real(r_double),intent(in) :: dt
      real(r_kind),intent(in):: zint(im,km+1),zdot(im,km+1)
      real(r_kind),intent(inout):: zadv(im,km,nm)
      integer :: nvcn
      real(r_kind) :: xvcn
      integer i,j,k,kk,n,ivcn(im)
      logical lvcn(im)
      real(r_kind) zdm,zda,zdb,vcn(im,km-1)
      real(r_kind) rnu,cm(im,km),cu(im,km-1),cl(im,km-1)
      real(r_kind) rr(im,km,nm)
      real(r_kind) cfl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute vertical courant number
!  increase by 10% for safety
! note: zdot goes top to bottom, other inputs go bottom to top
      cfl = 1.8
      nvcn=0
      xvcn=0.
      lvcn=.false.
      do k=1,km-1
        kk = km-k+2
        do i=1,im
          zdm=0.5*(zint(i,kk)-zint(i,kk+2))
          vcn(i,k)=abs(zdot(i,km-(k+1)+2)*dt/zdm)*cfl
          lvcn(i)=lvcn(i).or.vcn(i,k).gt.1
          xvcn=max(xvcn,vcn(i,k))
        enddo
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  determine points requiring filtering
      if(xvcn.gt.1) then
        do i=1,im
          if(lvcn(i)) then
            ivcn(nvcn+1)=i
            nvcn=nvcn+1
          endif
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute tridiagonal matrim
        do j=1,nvcn
          cm(j,1)=1
        enddo
        do k=1,km-1
          kk = km-k+2
          do j=1,nvcn
            i=ivcn(j)
            if(vcn(i,k).gt.1) then
              !zdm=zmid(i,k)-zmid(i,k+1)
              zdm=0.5*(zint(i,kk)-zint(i,kk+2))
              zda=zint(i,kk+1)-zint(i,kk+2)
              zdb=zint(i,kk)-zint(i,kk+1)
              rnu=(vcn(i,k)**2-1)/4
              cu(j,k)=-rnu*zdm/zdb
              cl(j,k)=-rnu*zdm/zda
              cm(j,k)=cm(j,k)-cu(j,k)
              cm(j,k+1)=1-cl(j,k)
            else
              cu(j,k)=0
              cl(j,k)=0
              cm(j,k+1)=1
            endif
          enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  fill fields to be filtered
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              rr(j,k,n)=zadv(i,k,n)
            enddo
          enddo
        enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  solve tridiagonal system
        call tridim_hyb(nvcn,im,km,km,nm,cl,cm,cu,rr,cu,rr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  replace filtered fields
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              zadv(i,k,n)=rr(j,k,n)
            enddo
          enddo
        enddo
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end

      subroutine tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!                .      .    .                                       .
! subprogram:    tridim_hyb      solves tridiagonal matrix problems.
!   prgmmr: iredell          org: w/nmc23    date: 91-05-07
!
! abstract: this routine solves multiple tridiagonal matrix problems
!   with multiple right-hand-side and solution vectors for every matrix.
!   the solutions are found by eliminating off-diagonal coefficients,
!   marching first foreward then backward along the matrix diagonal.
!   the computations are vectorized around the number of matrices.
!   no checks are made for zeroes on the diagonal or singularity.
!
! program history log:
!   97-07-30  iredell
!
! usage:    call tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!
!   input argument list:
!     l        - integer number of tridiagonal matrices
!     lx       - integer first dimension (lx>=l)
!     n        - integer order of the matrices
!     nx       - integer second dimension (nx>=n)
!     m        - integer number of vectors for every matrix
!     cl       - real(r_kind) (lx,2:n) lower diagonal matrix elements
!     cm       - real(r_kind) (lx,n) main diagonal matrix elements
!     cu       - real(r_kind) (lx,n-1) upper diagonal matrix elements
!                (may be equivalent to au if no longer needed)
!     r        - real(r_kind) (lx,nx,m) right-hand-side vector elements
!                (may be equivalent to a if no longer needed)
!
!   output argument list:
!     au       - real(r_kind) (lx,n-1) work array
!     a        - real(r_kind) (lx,nx,m) solution vector elements
!
! attributes:
!   language: fortran 77.
!   machine:  cray.
!
      use kinds, only: r_kind
      implicit none
      real(r_kind) cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m), &
                               au(lx,n-1),a(lx,nx,m)
      integer i,j,k,l,m,n,nx,lx
      real(r_kind) fk
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  march up
      do i=1,l
        fk=1./cm(i,1)
        au(i,1)=fk*cu(i,1)
      enddo
      do j=1,m
        do i=1,l
          fk=1./cm(i,1)
          a(i,1,j)=fk*r(i,1,j)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)=fk*cu(i,k)
        enddo
        do j=1,m
          do i=1,l
            fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
            a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
          enddo
        enddo
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  march down
      do j=1,m
        do i=1,l
          fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
          a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
        enddo
      enddo
      do k=n-1,1,-1
        do j=1,m
          do i=1,l
            a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
          enddo
        enddo
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
