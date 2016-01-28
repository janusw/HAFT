!*******************************************************************************
!  Read and interpolate Hades Acceptance Matrix for Theorists
!
!  Code version 2.0 of February 15, 2011
!
!  Created :  07/04/05 R. Holzmann
!  Modified : 25/04/05 R. Holzmann  added new interpolation modes
!  Modified : 04/05/05 R. Holzmann  added resolution function
!  Modified : 24/06/05 R. Holzmann  added pair acceptance
!  Modified : 20/07/05 R. Holzmann  fixed error messages on input
!  Modified : 06/02/06 R. Holzmann  support for diff. matrix sizes (v1.0) 
!  Modified : 15/04/08 R. Holzmann  added function to set resolution (v1.1) 
!  Modified : 07/02/09 R. Holzmann  added p-dependant Eloss in smear4momentum()
!  Modified : 16/04/09 R. Holzmann  theta and phi resolution from Ar+KCl embedding
!  Modified:  15/02/11 R. Holzmann  added support for non-gaussian momentum smearing 
!  Mofified:  20/04/11 J. Weil      use stream I/O, generic intrinsics, etc
!
!  Usage : 1) set acceptance file name with
!                 call setFileName(fname)
!          2) sample single acceptance values with calls to
!                 acc = getHadesAcceptance(id,mom,theta,phi,mode)
!          3) sample pair acceptance values with calls to
!                 acc = getHadesPairAcceptance(mass,pt,rapidity,mode)
!*******************************************************************************

module HAFT

  implicit none
  private


  public :: getHadesAcceptance, getHadesPairAcceptance
  public :: setFileName, setPairFileName
  public :: smearhadesmomentum, smearHades3Momentum, smearHadesPair
  public :: setResolutionParameters, setAngularResolutionParameters


  !  HAFT declaration of acceptance matrix arrays and resolution tables
  !  The dimensions MUST match all array sizes in the file!
  integer*4, parameter :: size  = 250000, &  ! <<== change if < xdim*ydim*zdim
                          sizep = 1000,   &  ! <<== change if < xdimp*ydimp
                          nids  = 14         ! <<== change if < max id

  real*4, parameter :: pi = 3.141592654, twopi = 2.*pi

  character(len=200), save :: fname  = 'HadesAcceptanceFilter.acc', &
                              fname2 = 'HadesPairAcceptanceFilter.acc'

  integer*4, save :: iseed = 123
  integer*4, save :: readflag = 0, readflag2 = 0

  character(len=80) :: comment, comment2
  integer*4, dimension(nids) :: xdim, ydim, zdim, resflag, xdimp, ydimp
  real*4, dimension(nids) :: dp, dth, dph, pmin, pmax, thmin, thmax, phmin, phmax
  real*4, dimension(nids) :: dpp, dthp, pminp, pmaxp, thminp, thmaxp

  ! matrices are declared for e+, e-, pi+, pi-, K+, K- and p
  real*4, dimension(size) :: matrix2, matrix3, matrix8, matrix9, matrix10, &
                             matrix12, matrix14, matrix51
  real*4 sigpA(3), sigpB(3), sigth, sigph, XX0
  integer*4 xdim2, ydim2, zdim2, ntab
  real*4 dm, dpt, drap, mmin, mmax, ptmin, ptmax, rapmin, rapmax

  real*4, dimension(sizep) :: par2p1, par2p2, par2p3, par2p4, par2p5, par2p6
  real*4, dimension(sizep) :: par3p1, par3p2, par3p3, par3p4, par3p5, par3p6


contains


  real*4 function getHadesAcceptance(pid,p0,theta0,phi0,mode)
    !
    !  Returns HADES acceptance for particle id of given momentum (in GeV/c),
    !  polar angle theta (in deg.) and azimuthal angle phi (in deg.)
    !  by table interpolation.
    !
    !  Lab frame used: z = beam axis, y = vertical axis
    !
    !                 ^ y
    !            x \  |
    !               \ |
    !                \|    z
    !                 O---->
    !
    !  id = 2 : positron
    !     = 3 : electron
    !     = 8 : pi+
    !     = 9 : pi-
    !     = 10: K+
    !     = 12: K-
    !     = 14: proton
    !
    !  mode = 0 : nearest-neighbour interpolation
    !       = 1 : tri-linear interpolation
    !       = 2 : tri-quadratic interpolation
    !       =-2 : tri-quadratic B-spline interpolation
    !       = 3 : tri-cubic interpolation
    !       =-3 : tri-cubic B-spline interpolation
    !       = 4 : tri-cubic Catmull-Rom spline
    !       =-4 : tri-cubic optimal cardinal spline
    !
      integer*4, intent(in) :: pid, mode
      real*4, intent(in) :: p0, theta0, phi0

      real*4 p, theta, phi
      real*4 u, v, w
      integer*4 ix, iy, iz, i, j, k
      integer*4 ilo, ihi, jlo, jhi, klo, khi
      real*4 sum, Kx, Ky, Kz, kernel
      integer*4 xdim, ydim, zdim
      real*4 plo, pup, dp, thlo, thup, dth, phlo, phup, dph

      integer*4 retcode, readHAFTmatrix
      real*4 getMatrixVal
      integer*4 mod_

      mod_ = -2   ! use B-spline interpolation
!     mod = mode

      getHadesAcceptance = 0.

      retcode = readHAFTmatrix()
      if (retcode.eq.-1) return

      call getDimensions(pid,xdim,ydim,zdim)
      call getLimits(pid,plo,pup,dp,thlo,thup,dth,phlo,phup,dph) 

      if (p0.lt.plo .or. theta0.lt.thlo .or. theta0.gt.thup) return
      if (phi0.lt.phlo) return

      p = min(p0,pup-2.01*dp)  ! level off acceptance at high p
      theta = theta0
      phi = phi0
      if (phi .gt. 60.) phi = mod(phi,60._4)   ! modulo HADES sector

      ix = xdim*((p-0.5*dp-plo)/(pup-plo)) + 1      ! floor indexes
      iy = ydim*((theta-0.5*dth-thlo)/(thup-thlo)) + 1
      iz = zdim*((phi-0.5*dph-phlo)/(phup-phlo)) + 1

      if (mod_.eq.0 .or. mod_.eq.1) then   ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      else if (abs(mod_).eq.2 .or. abs(mod_).eq.3 &
                              .or. abs(mod_).eq.4) then
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      else  ! mode not defined
        return
      end if

      if (ilo.lt.0 .or. jlo.lt.0 .or. klo.lt.0) return
      if (ihi.gt.xdim+1 .or. jhi.gt.ydim+1 .or. khi.gt.zdim+1) return

      sum = 0. 
      do i=ilo,ihi                      ! triple interpolation loop
        u = (p - (real(i)-0.5)*dp-plo)/dp
        Kx = kernel(u,mod_)
        do j=jlo,jhi
          v = (theta - (real(j)-0.5)*dth-thlo)/dth
          Ky = kernel(v,mod_)
          do k=klo,khi
            w = (phi - (real(k)-0.5)*dph-phlo)/dph
            Kz = kernel(w,mod_)
            sum = sum + getMatrixVal(i,j,k,pid)*Kx*Ky*Kz
          end do
        end do
      end do
      sum = max(min(sum, 1.01),0.)  ! clip over/undershoots

      getHadesAcceptance = sum

  end function getHadesAcceptance



  real*4 function getHadesPairAcceptance(mass0,pt0,rap0,mode)
    !
    !  Returns HADES pair acceptance for given mass (in GeV/c**2),
    !  transverse momentum (in GeV/c) and rapidity (in lab frame)
    !  by table interpolation.
    !
    !  Lab frame used: z = beam axis, y = vertical axis
    !
    !                 ^ y
    !            x \  |
    !               \ |
    !                \|    z
    !                 O---->
    !
    !
    !  mode = 0 : nearest-neighbour interpolation
    !       = 1 : tri-linear interpolation
    !       = 2 : tri-quadratic interpolation
    !       =-2 : tri-quadratic B-spline interpolation
    !       = 3 : tri-cubic interpolation
    !       =-3 : tri-cubic B-spline interpolation
    !       = 4 : tri-cubic Catmull-Rom spline
    !       =-4 : tri-cubic optimal cardinal spline
    !
      integer*4, intent(in) :: mode
      real*4, intent(in) :: mass0, pt0, rap0

      real*4 mass, pt, rap
      real*4 u, v, w
      integer*4 ix, iy, iz, i, j, k
      integer*4 ilo, ihi, jlo, jhi, klo, khi
      real*4 sum, Kx, Ky, Kz, kernel
      integer*4 xdim, ydim, zdim
      real*4 mlo, mup, dm, ptlo, ptup, dpt, raplo, rapup, drap

      logical executed
      data executed /.false./
      save executed
      integer*4 retcode, readHAFTPairmatrix
      real*4 getMatrixVal
      integer*4 mod

      getHadesPairAcceptance = 0.
      retcode = readHAFTPairMatrix()
      if (retcode.eq.-1) return

      mod = 1    ! use trilinear interpolation
!     mod = mode ! (use mode = 0 or 1, otherwise problems at pt=0!) 

      call getDimensions(51,xdim,ydim,zdim)
      call getLimits(51,mlo,mup,dm,ptlo,ptup,dpt,raplo,rapup,drap) 

      if (mass0.lt.mlo .or. pt0.lt.ptlo .or. pt0.gt.ptup &
          .or. rap0.lt.raplo .or. rap0.gt.rapup) return

      mass = min(mass0,mup-2.01*dm)  ! level off acceptance at high mass
      pt = pt0
      rap = rap0

      ix = xdim*((mass-0.5*dm-mlo)/(mup-mlo)) + 1      ! floor indexes
      iy = ydim*((pt-0.5*dpt-ptlo)/(ptup-ptlo)) + 1
      iz = zdim*((rap-0.5*drap-raplo)/(rapup-raplo)) + 1

      if (mod.eq.0 .or. mod.eq.1) then   ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      else if (abs(mod).eq.2 .or. abs(mod).eq.3 &
                             .or. abs(mod).eq.4) then
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      else  ! mode not defined
        return
      end if

      if (ilo.lt.0 .or. jlo.lt.0 .or. klo.lt.0) return
      if (ihi.gt.xdim+1 .or. jhi.gt.ydim+1 .or. khi.gt.zdim+1) return

      sum = 0. 
      do i=ilo,ihi                      ! triple interpolation loop
        u = (mass - (real(i)-0.5)*dm-mlo)/dm
        Kx = kernel(u,mod)
        do j=jlo,jhi
          v = (pt - (real(j)-0.5)*dpt-ptlo)/dpt
          Ky = kernel(v,mod)
          do k=klo,khi
            w = (rap - (real(k)-0.5)*drap-raplo)/drap
            Kz = kernel(w,mod)
            sum = sum + getMatrixVal(i,j,k,51)*Kx*Ky*Kz
          end do
        end do
      end do
      sum = max(min(sum, 1.01),0.)  ! clip over/undershoots

      getHadesPairAcceptance = sum

  end function getHadesPairAcceptance



  real*4 function kernel(u,mode)
    !
    !  Compute interpolation kernel
    !
    !  mode = 0: nearest neighbour
    !       = 1: piece-wise linear
    !       = 2: piece-wise quadratic
    !       =-2: quadratic B-spline
    !       = 3: cubic spline (B=0,C=1)
    !       =-3: cubic B-spline (B=1, C=0)
    !       = 4: cubic Catmull-Rom spline (B=0, C=1/2)
    !       =-4: cubic "optimal" cardinal spline (B=1/3, C=1/3)
    !
    !  mode >=0: interpolating (i.e. exact at grid points)
    !  mode < 0: approximating (i.e. not exact, but smoother)
    !
      real*4, intent(in) :: u
      integer*4, intent(in) :: mode

      real*4 ua

      select case (mode)
      case (0)  ! nearest neighbour
        if (u.gt.-0.5 .and. u.le.0.5) then
          kernel = 1.
        else 
          kernel = 0.
        end if

      case (1)  ! linear
        ua = abs(u)
        if (ua.lt.1.) then
          kernel = 1.-ua
        else 
          kernel = 0.
        end if

      case (2)  ! quadratic
        ua = abs(u)
        if (ua.le.0.5) then
          kernel = 1.0-2.*ua*ua
        else if (ua.le.1.5) then
          kernel = 1.5+(ua-2.5)*ua
        else 
          kernel = 0.
        end if

      case (-2)  ! quadratic B-spline
        ua = abs(u)
        if (ua.le.0.5) then
          kernel = 0.75-ua*ua
        else if (ua.le.1.5) then
          kernel = 1.125+0.5*(ua-3.)*ua
        else 
          kernel = 0.
        end if

      case (3)  ! cubic
        ua = abs(u)
        if (ua.le.1.) then
          kernel = 1.+(ua-2.)*ua*ua
        else if (ua.le.2.) then
          kernel = 4.+((5.-ua)*ua-8.)*ua
        else 
          kernel = 0.
        end if

      case (-3)  ! cubic B-spline
        ua = abs(u)
        if (ua.le.1.) then
          kernel = 2./3.+(0.5*ua-1.)*ua*ua
        else if (ua.le.2.) then
          kernel = 4./3.+((1.-ua/6.)*ua-2.)*ua
        else 
          kernel = 0.
        end if

      case (4)  ! cubic Catmull-Rom
        ua = abs(u)
        if (ua.le.1.) then
          kernel = 1.+(1.5*ua-2.5)*ua*ua
        else if (ua.le.2.) then
          kernel = 2.+((2.5-0.5*ua)*ua-4.)*ua
        else 
          kernel = 0.
        end if

      case (-4)  ! optimal cubic cardinal spline
                 ! (=compromise between blurring and ringing)
        ua = abs(u)
        if (ua.le.1.) then
          kernel = 8./9.+(7./6.*ua-2.)*ua*ua
        else if (ua.le.2.) then
          kernel = 16./9.+((2.-7./18.*ua)*ua-10./3.)*ua
        else 
          kernel = 0.
        end if

      case default  ! undefined mode
        kernel = 0.

      end select

  end function kernel



  integer*4 function readHAFTmatrix()
    !
    !  Opens file in unformatted direct access mode
    !  and reads HADES acceptance matrices (as linearized arrays)
    !
      integer*4, parameter :: runit = 77  ! change if input unit is already busy

      integer*4 pid
      integer*4 i
      integer*4 bytes       ! byte counter
      integer*4 bins
      integer*4 lc          ! string length

      readHAFTmatrix = 0

      if (readflag.eq.1) return

      readflag = 0
      do i=1,size           ! set all matrices to 0
        matrix2(i) = 0.
        matrix3(i) = 0.
        matrix8(i) = 0.
        matrix9(i) = 0. 
        matrix10(i) = 0.
        matrix12(i) = 0.
        matrix14(i) = 0.
      end do
      do i=1,nids
        xdim(i) = 0
        ydim(i) = 0
        zdim(i) = 0
        pmin(i) = 0.
        pmax(i) = 0.
        dp(i) = 0.
        thmin(i) = 0.
        thmax(i) = 0.
        dth(i) = 0.
        phmin(i) = 0.
        phmax(i) = 0.
        dph(i) = 0.
      end do

      lc = len(fname)
      do i=1,lc
         if (fname(i:i+3).eq.'.acc') goto 1
      end do
 1    lc = i+3
      open(unit=runit,file=fname,access='stream',status='old',err=99)
      bytes=1
      write(6,*) ' '
      read(runit,pos=bytes,err=100) comment
      write(6,'(a80)') comment
      write(6,*) '--------------------------------------'
      bytes = bytes + 80
      read(runit,pos=bytes,err=100) sigpA(1), sigpA(2), sigpA(3), &
                                    sigpB(1), sigpB(2), sigpB(3), &
                                    sigth, sigph, XX0
      bytes = bytes + 9*4

 40   read(runit,pos=bytes,end=50,err=100) pid  ! break out if EOF reached
      if (pid.lt.0) goto 45
      if (pid.lt.1 .or. pid.gt.nids) goto 50

!ccccccccccccccc read acceptance matrix 

      bytes = bytes + 4
      read(runit,pos=bytes,err=100) xdim(pid), ydim(pid), zdim(pid)

      bins = xdim(pid)*ydim(pid)*zdim(pid)
      if (bins.gt.size) goto 101  ! check if enough memory available

      bytes = bytes + 3*4
      read(runit,pos=bytes,err=100) pmin(pid),pmax(pid),   &
                                    thmin(pid),thmax(pid), &
                                    phmin(pid),phmax(pid)
      bytes = bytes + 6*4
      if (pid.eq.2) then        ! positron
        read(runit,pos=bytes,err=100) (matrix2(i),i=1,bins)
        write(6,'(''Matrix read for e+'')')
      else if (pid.eq.3) then     ! electron
        read(runit,pos=bytes,err=100) (matrix3(i),i=1,bins)
        write(6,'(''Matrix read for e-'')')
      else if (pid.eq.8) then     ! pi+
        read(runit,pos=bytes,err=100) (matrix8(i),i=1,bins)
        write(6,'(''Matrix read for pi+'')')
      else if (pid.eq.9) then     ! pi-
        read(runit,pos=bytes,err=100) (matrix9(i),i=1,bins)
        write(6,'(''Matrix read for pi-'')')
      else if (pid.eq.10) then    ! K+
        read(runit,pos=bytes,err=100) (matrix10(i),i=1,bins)
        write(6,'(''Matrix read for K+'')')
      else if (pid.eq.12) then    ! K-
        read(runit,pos=bytes,err=100) (matrix12(i),i=1,bins)
        write(6,'(''Matrix read for K-'')')
      else if (pid.eq.14) then    ! proton
        read(runit,pos=bytes,err=100) (matrix14(i),i=1,bins)
        write(6,'(''Matrix read for p'')')
      else
        write(6,'(''Unsupported particle ID: '',I3,''  STOP!'')')
        stop
      end if
      resflag(pid) = 0
      bytes = bytes + bins*4

!     write(6,*) 'Acceptance matrix for ID= ',pid
      write(6,*) 'dims= ',xdim(pid), ' ', ydim(pid), ' ', zdim(pid)
      write(6,*) 'lims= ',pmin(pid), ' ', pmax(pid), ' ', thmin(pid), &
                 ' ', thmax(pid), ' ', phmin(pid), ' ', phmax(pid)
      write(6,*) 'size of matrix :', bins
      write(6,*) '--------------------------------------'
      dp(pid) = (pmax(pid)-pmin(pid))/real(xdim(pid))
      dth(pid) = (thmax(pid)-thmin(pid))/real(ydim(pid))
      dph(pid) = (phmax(pid)-phmin(pid))/real(zdim(pid))

      goto 40  ! loop until EOF is reached

!ccccccccccccccc read resolution parameters 

 45   pid = -pid
      if (pid.lt.2 .or. pid.gt.3) goto 102

      bytes = bytes + 4
      read(runit,pos=bytes,err=100) xdimp(pid), ydimp(pid)

      bins = xdimp(pid)*ydimp(pid)
      if (bins.gt.sizep) goto 101  ! check if enough memory available

      bytes = bytes + 2*4
      read(runit,pos=bytes,err=100) pminp(pid),pmaxp(pid), &
                                    thminp(pid),thmaxp(pid)
      bytes = bytes + 4*4
      read(runit,pos=bytes,err=100) ntab ! nb. of parameter tables
      bytes = bytes + 4
      if (pid.eq.2) then          ! positron
        read(runit,pos=bytes,err=100) (par2p1(i),i=1,bins)
        bytes = bytes + bins*4
        if (ntab.gt.1) then
          read(runit,pos=bytes,err=100) (par2p2(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.2) then
          read(runit,pos=bytes,err=100) (par2p3(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.3) then
          read(runit,pos=bytes,err=100) (par2p4(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.4) then
          read(runit,pos=bytes,err=100) (par2p5(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.5) then
          read(runit,pos=bytes,err=100) (par2p6(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        write(6,'(''Resolution tables read for e+'')')
      else if (pid.eq.3) then     ! electron
        read(runit,pos=bytes,err=100) (par3p1(i),i=1,bins)
        bytes = bytes + bins*4
        if (ntab.gt.1) then
          read(runit,pos=bytes,err=100) (par3p2(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.2) then
          read(runit,pos=bytes,err=100) (par3p3(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.3) then
          read(runit,pos=bytes,err=100) (par3p4(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.4) then
          read(runit,pos=bytes,err=100) (par3p5(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab.gt.5) then
          read(runit,pos=bytes,err=100) (par3p6(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        write(6,'(''Resolution tables read for e-'')')
      else
        write(6,'(''Unsupported PID: '',I3,'' Use default smearing!'')')
        resflag(pid) = 0
        goto 40
      end if
      resflag(pid) = 1

!     write(6,*) 'Parameter tables for ID= ',pid
      write(6,*) 'dims= ',xdimp(pid), ' ', ydimp(pid)
      write(6,*) 'lims= ',pminp(pid), ' ', pmaxp(pid), ' ', thminp(pid), &
                 ' ', thmaxp(pid)
      write(6,*) 'size of parameter tables :', ntab, ' x', bins
      write(6,*) '--------------------------------------'
      dpp(pid) = (pmaxp(pid)-pminp(pid))/real(xdimp(pid))
      dthp(pid) = (thmaxp(pid)-thminp(pid))/real(ydimp(pid))

      goto 40  ! loop until EOF is reached

 50   close(runit)


      readHAFTmatrix = bytes-1 ! return number of bytes read
      readflag = 1
      return

!     Error opening or reading

 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',fname(1:lc)
      readHAFTMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',fname(1:lc)
      readHAFTmatrix = -1
      return
 101  close(runit)
      write(6,*) 'Size error: ', bins, ' >', size,' File = ',fname(1:lc)
      readHAFTmatrix = -1
      return
 102  close(runit)
      write(6,*) 'PID not yet supported: ', pid, ' File = ',fname(1:lc)
      readHAFTmatrix = -1
      return
  end function readHAFTmatrix



  integer*4 function readHAFTPairMatrix()
    !
    !  Opens file in unformatted direct access mode
    !  and reads HADES pair acceptance matrix (as linearized array)
    !
      integer*4, parameter :: runit = 78  ! change if input unit is already busy

      integer*4 i
      integer*4 bytes       ! byte counter
      integer*4 bins
      integer*4 lc

      readHAFTPairMatrix = 0

      if (readflag2.eq.1) return

      readflag2 = 0
      do i=1,size           ! set matrix to 0
        matrix51(i) = 0.
      end do

      lc = len(fname2)
      do i=1,lc
         if (fname2(i:i+3).eq.'.acc') goto 1
      end do
 1    lc = i+3
      open(unit=runit,file=fname2,access='stream',status='old',err=99)
      bytes=1
      read(runit,pos=bytes,err=100) comment2
      write(6,'(a80)') comment2
      bytes = bytes + 80
      read(runit,pos=bytes,err=100) xdim2, ydim2, zdim2

      bins = xdim2*ydim2*zdim2
      if (bins.gt.size) goto 101  ! check if enough memory available

      bytes = bytes + 3*4
      read(runit,pos=bytes,err=100) mmin,mmax,ptmin,ptmax,rapmin,rapmax
      bytes = bytes + 6*4
      read(runit,pos=bytes,err=100) (matrix51(i),i=1,bins)
      write(6,'(''Matrix read for e+e- pairs'')')
      bytes = bytes + bins*4
      close(runit)

      dm = (mmax-mmin)/real(xdim2)
      dpt = (ptmax-ptmin)/real(ydim2)
      drap = (rapmax-rapmin)/real(zdim2)

      write(6,'('' coms= '',a80)') comment2
      write(6,*) 'dims= ',xdim2, ' ', ydim2, ' ', zdim2
      write(6,*) 'lims= ',mmin, ' ', mmax, ' ', ptmin, ' ', ptmax, &
                 ' ', rapmin, ' ', rapmax
      write(6,*) 'size of matrix :', bins

      readHAFTPairMatrix = bytes-1 ! return number of bytes read
      readflag2 = 1
      return

!     Error opening or reading

 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',fname2(1:lc)
      readHAFTPairMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',fname2(1:lc)
      readHAFTPairMatrix = -1
      return
 101  close(runit)
      write(6,*) 'Size error: ', bins,' >',size, ' File = ',fname2(1:lc)
      readHAFTPairMatrix = -1
      return
  end function readHAFTPairMatrix



  subroutine setFileName(name)
    !
    !  Sets name of input file containing the filter
    !
      character*(*), intent(in) :: name

      integer*4 dummy, readHAFTmatrix

      fname = name
      dummy = readHAFTmatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname |'',a80,''|'')') fname

  end subroutine setFileName



  subroutine setPairFileName(name)
    !
    !  Sets name of input file containing the filter
    !
      character*(*), intent(in) :: name

      integer*4 dummy, readHAFTPairMatrix

      fname2 = name
      dummy = readHAFTPairMatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname2 |'',a80,''|'')') fname2

  end subroutine setPairFileName



  real*4 function getMatrixVal(i,j,k,pid)
    !
    !  Returns acceptance value at cell (i,j,k) of linearized matrix
    !  for particle ID
    !
      integer*4, intent(in) :: i, j, k, pid

      integer*4 xdi, ydi, zdi
      integer*4 i1, j1, k1, ilin

      call getDimensions(pid,xdi,ydi,zdi)

      i1 = min(max(1,i),xdi)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),ydi)  ! This effectively extrapolates matrix
      k1 = min(max(1,k),zdi)  ! beyond table boundaries.

      ilin = i1+xdi*(j1-1)+xdi*ydi*(k1-1)  ! linearized index
      select case (pid)
      case (2)     ! positron
        getMatrixVal = matrix2(ilin)
      case (3)     ! electron
        getMatrixVal = matrix3(ilin)
      case (8)     ! pi+
        getMatrixVal = matrix8(ilin)
      case (9)     ! pi-
        getMatrixVal = matrix9(ilin)
      case (10)    ! K+
        getMatrixVal = matrix10(ilin)
      case (12)    ! K-
        getMatrixVal = matrix12(ilin)
      case (14)    ! proton
        getMatrixVal = matrix14(ilin)
      case (51)    ! dilepton
        getMatrixVal = matrix51(ilin)
      case default
        getMatrixVal = 0.
      end select
  end function getMatrixVal



  subroutine getLimits(pid,xlo,xhi,dx,ylo,yhi,dy,zlo,zhi,dz) 
    !
    !  Return the lower and upper limits, and step sizes of the table
    !
      integer*4, intent(in) :: pid
      real*4, intent(out) :: xlo, xhi, dx, ylo, yhi, dy, zlo, zhi, dz

      if (pid.eq.51) then  ! pair
        xlo  = mmin
        xhi  = mmax
        dx   = dm
        ylo = ptmin
        yhi = ptmax
        dy  = dpt
        zlo = rapmin
        zhi = rapmax
        dz = drap
      else
        xlo  = pmin(pid)
        xhi  = pmax(pid)
        dx   = dp(pid)
        ylo = thmin(pid)
        yhi = thmax(pid)
        dy  = dth(pid)
        zlo = phmin(pid)
        zhi = phmax(pid)
        dz = dph(pid)
      end if
  end subroutine getLimits



  subroutine getDimensions(pid,nx,ny,nz)
    !
    !  Return the dimensions of a table of particle pid
    !
      integer*4, intent(in)  :: pid
      integer*4, intent(out) :: nx, ny, nz

      if (pid.eq.51) then  ! pair
        nx  = xdim2
        ny = ydim2
        nz = zdim2
      else
        nx = xdim(pid)
        ny = ydim(pid)
        nz = zdim(pid)
      end if
  end subroutine getDimensions



  subroutine smearHades4Momentum(mom,mode,pid)
    !
    !  Apply the Hades momentum resolution to a 4-momentum vector (in GeV/c)
    !
    !  Lab frame used: z = beam axis, y = vertical axis
    !
    !                 ^ y
    !            x \  |
    !               \ |
    !                \|    z
    !                 O---->
    !
    !  All components of the 4-momentum vector are changed by this call.
    !
    !  If parameter tables are loaded, the particle momentum is smeared according
    !  to an asymmetric response function, if not, default gaussian sampling is used.
    !
    !  The default resolution mode is determined by:
    !
    !   mode = 1 : low-resolution     (MDC 1+2)
    !        = 2 : medium-resolution  (MDC 1+2+3)
    !        = 3 : high-resolution    (MDC 1+2+3+4)
    !
    !
      real*4, intent(inout) :: mom(4)
      integer*4, intent(in) :: mode,pid

      integer*4 retcode, readHAFTmatrix
      real*4 mass, mass2, pt, pt2, ptot, ptot2, theta, phi, sinth
      real*4 sigp, sampleGauss, betainv, sigms, sigms2, sigthms, sigphms
      real*4, parameter :: r2d = 180./pi
      real*4 ploss
      real*4 respar(10), param, sampleMP
      integer*4 i

      if (readflag.eq.0) then
        retcode = readHAFTmatrix()
        if (retcode.eq.-1) return
      end if

      pt2 = mom(1)**2 + mom(2)**2
      pt = sqrt(pt2)
      ptot2 = pt2 + mom(3)**2
      ptot = sqrt(ptot2)                      ! total momentum in GeV/c
      mass2 = mom(4)**2 - ptot2               ! particle mass does not change
      if (mass2.lt.0.) mass2 = 0.
      mass = sqrt(mass2)
      betainv = 1.
      if (mass2.gt.0.) betainv = sqrt(1. + mass2/ptot2)  ! 1/beta
      sigms = 0.0136*betainv/ptot*sqrt(XX0)*(1+0.038*log(XX0))
                                              ! multiple scattering angle
      sigms2 = sigms*sigms
      sigthms = sqrt(sigth*sigth+sigms2)      ! add quadratically to resolution
      sigphms = sqrt(sigph*sigph+sigms2)

      if (pid.eq.2 .or. pid.eq.3) then  ! leptons, use Ar+KCl embedding values
         sigthms = 0.00246/ptot + 0.00045
         sigphms = 0.00294/ptot + 0.00059
      endif

      theta = acos(mom(3)/ptot)               ! polar angle in radian
      sinth = sin(theta)
      phi = 0.
      if (pt.gt.0.) phi = acos(mom(1)/pt)
      if (mom(2).lt.0.) phi = twopi - phi     ! azimuthal angle in radian

!  If resolution parameters are available, use dedicated smearing

      if (resflag(pid).eq.1) then   ! resolution parameters are loaded for pid?

!        write(6,*) mom(1), mom(2), mom(3), mass, ptot, r2d*theta
         do i=1,ntab    ! look up parameters from tables
            respar(i) = param(ptot,r2d*theta,pid,i)
!           write(6,*) i, respar(i)
         end do

         ptot = ptot*(1.+sampleMP(respar,2._4))  ! randomize momentum

      else  ! use default gaussian smearing

         if (mode.lt.1 .or. mode.gt.3) return  ! unknown mode
         sigp = 0.01*ptot*(sigpA(mode)+sigpB(mode)*ptot)  ! momentum resolution
         ptot = max(0.,sampleGauss(ptot,sigp))  ! smear total momentum

         if (pid.eq.2 .or. pid.eq.3) then
            ploss = 0.0018
         else
            ploss = 0.0
         end if
         if (ptot.gt.0.1) then
            if (pid.eq.2) &
                  ploss = ploss + 0.005*(1.-exp(-(ptot-0.1)/0.053)) ! positron
            if (pid.eq.3) &
                  ploss = ploss + 0.006*(1.-exp(-(ptot-0.1)/0.078)) ! electron
         end if
         ptot = ptot - ploss                    ! Eloss of electrons in target

      end if

      theta = abs(sampleGauss(theta,sigthms)) ! smear polar angle
      if (theta.gt.pi) theta = twopi - theta  ! if > pi, mirror angle
      if (sinth.gt.0.01) phi = sampleGauss(phi,sigphms/sinth) ! smear azimuth
      if (phi.lt.0.) phi = phi + twopi        ! and check if within range 
      if (phi.gt.twopi) phi = phi - twopi

      sinth = sin(theta)
      mom(1) = ptot*sinth*cos(phi)            ! new momentum components
      mom(2) = ptot*sinth*sin(phi)
      mom(3) = ptot*cos(theta)
      mom(4) = sqrt(mom(1)**2 + mom(2)**2 + mom(3)**2 + mass2)  ! total energy

  end subroutine smearHades4Momentum



  subroutine smearHadesMomentum(p,mode,pid)
      real*8, intent(inout) :: p(0:3)
      integer*4, intent(in) :: mode, pid

      real*4 mom4(4)

      mom4(1:3) = p(1:3)
      mom4(4) = p(0)

      call smearHades4Momentum(mom4,mode,pid)

      p(1:3) = mom4(1:3)
      p(0) = mom4(4)
  end subroutine smearHadesMomentum



  subroutine smearHades3Momentum(mom3,mode,pid)
    !
    !  Apply Hades momentum resolution to a 3-momentum (calculate multiple
    !  scattering assuming the particle is an electron) 
    !
      real*4, intent(inout) :: mom3(3)
      integer*4, intent(in) :: mode,pid

      real*4 mom4(4)
      real*4 mass

      if (mode.lt.1 .or. mode.gt.3) return    ! unknown mode

      if (pid.eq.2 .or. pid.eq.3) then        ! e+ or e-
         mass = 0.000510998918
      else if (pid.eq.8 .or. pid.eq.9) then   ! pi+ or pi- 
         mass = 0.13957018
      else if (pid.eq.11 .or. pid.eq.12) then ! K+ or K- 
         mass = 0.493677
      else if (pid.eq.14) then                ! proton
         mass = 0.938272029
      else
         return  ! particle is not supported
      end if

      mom4(1) = mom3(1)
      mom4(2) = mom3(2)
      mom4(3) = mom3(3)
      mom4(4) = sqrt(mom3(1)**2+mom3(2)**2+mom3(3)**2 + mass**2)
      call smearHades4Momentum(mom4,mode,pid)
      mom3(1) = mom4(1)
      mom3(2) = mom4(2)
      mom3(3) = mom4(3)
  end subroutine smearHades3Momentum



  subroutine smearHadesPair(pair,mode)
    !
    !  Apply Hades momentum resolution to a pair (calculate multiple
    !  scattering assuming the particle is an electron) 
    !
      real*4, intent(inout) :: pair(3)
      integer*4, intent(in) :: mode

      real*4 m, pt, rap, sigpt, sigm, sigrap
      real*4 sampleGauss
      integer*4 retcode, readHAFTmatrix

      real*4 mtab(10), par1(10), par2(10), par3(10), par4(10), par5(10)
      data mtab /0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2/  ! mass grid
      data par1 /0.0077, -0.0082, -0.0125, -0.0120, -0.0114, &
                -0.0106, -0.0098, -0.0085, -0.0078, -0.0075/
      data par2 /0.0820, 0.0460, 0.0260, 0.0210, 0.0190, &
                 0.0183, 0.0182, 0.0181, 0.0180, 0.0180/
      data par3 /13.5, 16.2, 19.9, 20.2, 19.2, &
                 18.0, 16.9, 14.8, 12.8, 11.0/
      data par4 /-10.0, -15.4, -20.4, -23.1, -22.3, &
                 -21.6, -20.6, -19.4, -18.2, -17.9/
      data par5 /18.1, 11.6, 10.4, 10.0, 9.4, 8.5, 7.8, 7.0, 6.2, 5.7/

      real*4 respar(10), interpol, sampleMP

      if (readflag.eq.0) then
        retcode = readHAFTmatrix()
        if (retcode.eq.-1) return
      end if

      m = pair(1)
      pt = pair(2)
      rap = pair(3)

      if (mode.eq.4) then  ! use skewd mass function
        sigpt = 0.01*pt*(sigpA(3)+sigpB(3)*pt)/sqrt(2.)       ! pt resolution
        pt = max(0.,sampleGauss(pt,sigpt))                  ! smear pt

        respar(1) = interpol(m,mtab,par1,10)
        respar(2) = interpol(m,mtab,par2,10)
        respar(3) = interpol(m,mtab,par3,10)
        respar(4) = interpol(m,mtab,par4,10)
        respar(5) = interpol(m,mtab,par5,10)
!        do i=1,5    ! interpolate parameters
!          write(6,*) i, respar(i)
!        end do
        m = m*(1.+sampleMP(respar,sqrt(2._4)))                  ! randomize mass
      else ! use gaussian smearing
        sigpt = 0.01*pt*(sigpA(mode)+sigpB(mode)*pt)/sqrt(2.) ! pt resolution
        pt = max(0.,sampleGauss(pt,sigpt))                  ! smear pt

        sigm = 0.01*m*(sigpA(mode)+sigpB(mode)*m)/sqrt(2.)
        m = max(0.,sampleGauss(m,sigm))                     ! smear mass
      end if

      sigrap = 0.1    !  just a guess!
      rap = sampleGauss(rap,sigrap)

      pair(1) = m
      pair(2) = pt
      pair(3) = rap
  end subroutine smearHadesPair



  real*4 function sampleGauss(mean,sigma)
    !
    !  Return random number according to a normal distribution.
    !
    !  Calls ran(iseed), a uniform random number generator returning ]0,1[.
    !
      real*4, intent(in) :: mean, sigma

      real*4 theta

      sampleGauss = mean
      if (sigma.le.0.) return
      theta = twopi*ran(iseed)
      sampleGauss = mean + sigma*cos(theta)*sqrt(-2.*log(ran(iseed)))
  end function sampleGauss



  real*4 function momSpread(x,respar,ns)
    !
    !     HADES momentum spread (pRec-pSim)/pSim
    !
      real*4, intent(in) :: x, respar(10), ns

      real*4 pos, sig, left, right, farleft
      real*4 argn, argp, argn2
      real*4 e2, lg10
      real*4 amp

      e2 = exp(-0.5*ns*ns)
      lg10 = log(10.)

      pos = respar(1)     ! Mean
      sig = respar(2)     ! Sigma
      left = respar(3)    ! Par3 (>0)
      right = respar(4)   ! Par4 (<0)
      farleft = respar(5) ! Par5 (>0)

      if (x.ge.(pos-ns*sig) .or. x.le.(-lg10/left+pos-ns*sig)) then
         argn = 0.
      else 
         argn = 1.
      end if

      if (x.ge.(pos-ns*sig)) then
        argp = 1.
      else
        argp = 0.
      end if

      if (x.gt.(-lg10/left+pos-ns*sig)) then
        argn2 = 0.
      else
        argn2 = 1.
      end if

      amp = e2  ! Gauss amplitude at +/-2 sigma

      momSpread = exp( -0.5*((x-pos)/sig)*((x-pos)/sig) )               &  ! Gauss
             + amp*exp(  left*(x-(pos-ns*sig)) )*argn                   &  ! left tail (connects to Gauss at pos-ns*sig
             + amp*exp( right*(x-(pos-ns*sig)) )*argp                   &  ! right tail (Gauss sits on top of it)
             + 0.1*amp*exp( farleft*(x-(-lg10/left+pos-ns*sig)) )*argn2 &  ! far left tail
                                                                        &  ! (joins left tail where decayed to 1/10)
  end function momSpread



  real*4 function sampleMP(respar,ns)
    !
    !  Return random number according to the normalized HADES momentum distribution.
    !
    !  Calls ran(iseed), a uniform random number generator returning ]0,1[.
    !
      real*4, intent(in) :: respar(10), ns

      real*4 pos, sig, left, right, farleft
      real*4 A0, A1, A2, A3
      real*4 F0, F1, F2, F3, F
      real*4 dx, ftest
      real*4 r1, r2, r3
      integer*4 cnt, cnt1, cnt2, cnt3
      real*4 momSpread
      real*4 e2, lg10

      e2 = exp(-0.5*ns*ns)
      lg10 = log(10.)

      pos = respar(1)      ! centroid
      sig = respar(2)      ! width
      left = respar(3)     ! left slope
      right = respar(4)    ! right slope
      farleft = respar(5)  ! far left slope

!    compute function amplitudes
      A0 = 0.1*(1. + e2)
      A1 = 1. + e2
      A2 = 1. + exp(right*ns*sig)
      A3 = (e2 + exp(right*2.*ns*sig))/exp(right*2.*ns*sig)

!    compute function areas
      F0 = A0/farleft * (1. - exp(farleft*(lg10/left+ns*sig-pos-1.)))  ! [-1, 1/10left]
      F1 = A1/left * 9./10.                                            ! [1/10, pos-ns*sig]
      F2 = A2 * 2.*ns*sig                                              ! [pos-ns*sig, pos+ns*sig]
      F3 = A3/right * (exp(right*(1.-pos+ns*sig))-exp(right*2.*ns*sig))! [pos + ns*sig, 1]

      F = F0 + F1 + F2 + F3
      F0 = F0/F   ! normalize areas
      F1 = F1/F
      F2 = F2/F
      F3 = F3/F

!    sample dx by comparing with piece-wise function

      do cnt=1,1000    ! allow max 1000 trials
         cnt1 = 0
         cnt2 = 0
         cnt3 = 0
         r1 = ran(iseed)
!    select region and sample test function
         if  (r1.lt.F0) then             ! far left tail

 10         continue
               cnt1 = cnt1 + 1
               r2 = log(ran(iseed))
               dx = r2/farleft - ns*sig + pos - lg10/left
            if (cnt1.eq.1000) write(6,*) 'cnt1=1000 ', pos, sig, farleft
            if (dx.lt.-1. .and. cnt1.lt.1000) goto 10  ! limit range to >=-1
            ftest = A0 * exp(farleft*(dx+lg10/left-pos+ns*sig))

         else if (r1.lt.F0+F1) then      ! left tail

 20         continue
               cnt2 = cnt2 + 1
               r2 = log(ran(iseed))
               dx = r2/left - ns*sig + pos
            if (cnt2.eq.1000) write(6,*) 'cnt2=1000', pos, sig, left
            if (dx.lt.-lg10/left-ns*sig+pos .and. cnt2.lt.1000) goto 20
            ftest = A1 * exp(left*(dx-pos+ns*sig))

         else if (r1.lt.F0+F1+F2) then   ! peak region

            r2 = ran(iseed) - 0.5
            dx = 2.*ns*sig*r2 + pos
            ftest = A2

         else                            ! right tail

 30         continue
               cnt3 = cnt3 + 1
               r2 = log(ran(iseed))
               dx = r2/right + ns*sig + pos
            if (cnt3.eq.1000) write(6,*) 'cnt3=1000', pos, sig, right
            if (dx.gt.1. .and. cnt3.lt.1000) goto 30 ! limit range to <=1
            ftest = A3 * exp(right*(dx-pos+ns*sig))

         end if

!   do rejection test
         sampleMP = dx

         r3 = ran(iseed)
         if ( r3.lt.momSpread(dx,respar,ns)/ftest ) return

      end do

      write(6,*) 'cnt=1000'
      sampleMP = 0.

  end function sampleMP



  subroutine setResolutionParameters(mode,a,b)
    !
    !     Set momentum resolution parameters
    !
      real*4, intent(in) :: a, b
      integer*4, intent(in) :: mode

      if (mode.ge.1 .and. mode.le.3) then
         sigpA(mode) = a
         sigpB(mode) = b
      endif

      write(6,*) 'mode=1   ', ' ', sigpA(1), ' ', sigpB(1)
      write(6,*) 'mode=2   ', ' ', sigpA(2), ' ', sigpB(2)
      write(6,*) 'mode=3   ', ' ', sigpA(3), ' ', sigpB(3)
  end subroutine setResolutionParameters



  subroutine setAngularResolutionParameters(stheta,sphi)
    !
    !     Set angular resolution parameters
    !
      real*4, intent(in) :: stheta, sphi

      if (stheta.gt.0. .and. sphi.gt.0.) then
         sigth = stheta
         sigph = sphi
      endif

      write(6,*) 'dTh, dPh:', ' ', sigth, ' ', sigph
  end subroutine setAngularResolutionParameters



  real*4 function param(pin,thin,pid,itab)
    !
    !     Interpolate resolution parameter table as function
    !     of momentum and theta (pin in GeV/c and theta in degree)
    !
      real*4, intent(in) :: pin, thin
      integer*4, intent(in) :: pid, itab

      real*4 p, th
      integer*4 xdi, ydi
      real*4 plo, pup, dp0, thlo, thup, dth0
      integer*4 i, j, ix, iy, ilo, ihi, jlo, jhi
      integer*4 mod
      real*4 sum, u, v, kernel, Kx, Ky, getTableVal

      mod = 1
      param = 0.
      if (pid.lt.1 .or. pid.gt.nids) return

      p = pin
      th = thin
      plo = pminp(pid)
      pup = pmaxp(pid)
      dp0 = dpp(pid)
      thlo = thminp(pid)
      thup = thmaxp(pid)
      dth0 = dthp(pid)
      if (p.lt.plo) p = plo   ! safety fence
      if (p.gt.pup) p = pup
      if (th.lt.thlo) th = thlo
      if (th.gt.thup) th = thup

      xdi = xdimp(pid) ! get table dimensions
      ydi = ydimp(pid)

      ix = xdi*((p-0.5*dp0-plo)/(pup-plo)) + 1      ! floor indices
      iy = ydi*((th-0.5*dth0-thlo)/(thup-thlo)) + 1

      if (mod.eq.0 .or. mod.eq.1) then   ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
      else if (abs(mod).eq.2 .or. abs(mod).eq.3 &
                             .or. abs(mod).eq.4) then
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
      else  ! mode not defined
        return
      end if

      if (ilo.lt.0 .or. jlo.lt.0) return
      if (ihi.gt.xdi+1 .or. jhi.gt.ydi+1) return

      sum = 0. 
      do i=ilo,ihi                      ! double interpolation loop
        u = (p - (real(i)-0.5)*dp0-plo)/dp0
        Kx = kernel(u,mod)
        do j=jlo,jhi
          v = (th - (real(j)-0.5)*dth0-thlo)/dth0
          Ky = kernel(v,mod)
          sum = sum + getTableVal(i,j,pid,itab)*Kx*Ky
        end do
      end do

      param = sum
  end function param



  real*4 function getTableVal(i,j,pid,itab)
    !
    !  Returns acceptance value at cell (i,j) of linearized
    !  parameter table for particle ID
    !
      integer*4, intent(in) :: i, j, pid, itab

      integer*4 xdi, ydi
      integer*4 i1, j1, ilin

      getTableVal = 0.
      if (pid.lt.1 .or. pid.gt.nids) return

      xdi = xdimp(pid) ! get table dimensions
      ydi = ydimp(pid)

      i1 = min(max(1,i),xdi)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),ydi)  ! This effectively extrapolates matrix

      ilin = i1+xdi*(j1-1)    ! linearized index

      if (pid.eq.2) then          ! positron
        if (itab.eq.1) getTableVal = par2p1(ilin)
        if (itab.eq.2) getTableVal = par2p2(ilin)
        if (itab.eq.3) getTableVal = par2p3(ilin)
        if (itab.eq.4) getTableVal = par2p4(ilin)
        if (itab.eq.5) getTableVal = par2p5(ilin)
        if (itab.eq.6) getTableVal = par2p6(ilin)
      else if (pid.eq.3) then     ! electron
        if (itab.eq.1) getTableVal = par3p1(ilin)
        if (itab.eq.2) getTableVal = par3p2(ilin)
        if (itab.eq.3) getTableVal = par3p3(ilin)
        if (itab.eq.4) getTableVal = par3p4(ilin)
        if (itab.eq.5) getTableVal = par3p5(ilin)
        if (itab.eq.6) getTableVal = par3p6(ilin)
      end if
  end function getTableVal



  real*4 function interpol(x,xtab,ytab,n)
    !
    !     linear interpolation in table (xtab,ytab) 
    !
      real*4, intent(in) :: x, xtab(*), ytab(*)
      integer*4, intent(in) :: n

      integer*4 i
      real*4 a, b

      if (x.le.xtab(1)) then ! below table range
        interpol = ytab(1)
        return
      else if (x.ge.xtab(n)) then ! above table range
        interpol = ytab(n)
        return
      end if

      do i=2,n
        interpol = ytab(i)
        if (x.eq.xtab(i)) return
        if (x.lt.xtab(i)) goto 10
      end do
 10   continue

      a = ytab(i-1)
      b = (ytab(i)-ytab(i-1))/(xtab(i)-xtab(i-1))
      interpol = a + (x-xtab(i-1))*b
  end function interpol


end module
