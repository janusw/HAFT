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
  integer(4), parameter :: size  = 250000, &  ! <<== change if < xdim*ydim*zdim
                          sizep = 1000,   &  ! <<== change if < xdimp*ydimp
                          nids  = 14         ! <<== change if < max id

  real(4), parameter :: pi = 3.141592654, twopi = 2.*pi, lg10 = log(10.)

  character(len=200), save :: fname  = 'HadesAcceptanceFilter.acc', &
                              fname2 = 'HadesPairAcceptanceFilter.acc'

  integer(4), save :: readflag = 0, readflag2 = 0

  character(len=80) :: comment, comment2
  integer(4), dimension(nids) :: xdim, ydim, zdim, resflag, xdimp, ydimp
  real(4), dimension(nids) :: dp, dth, dph, pmin, pmax, thmin, thmax, phmin, phmax
  real(4), dimension(nids) :: dpp, dthp, pminp, pmaxp, thminp, thmaxp

  ! matrices are declared for e+, e-, pi+, pi-, K+, K- and p
  real(4), dimension(size) :: matrix2, matrix3, matrix8, matrix9, matrix10, &
                             matrix12, matrix14, matrix51
  real(4) sigpA(3), sigpB(3), sigth, sigph, XX0
  integer(4) xdim2, ydim2, zdim2, ntab
  real(4) dm, dpt, drap, mmin, mmax, ptmin, ptmax, rapmin, rapmax

  real(4), dimension(sizep) :: par2p1, par2p2, par2p3, par2p4, par2p5, par2p6
  real(4), dimension(sizep) :: par3p1, par3p2, par3p3, par3p4, par3p5, par3p6


contains


  !*****************************************************************************
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
  !*****************************************************************************
  real(4) function getHadesAcceptance(pid,p0,theta0,phi0,mode)

      integer(4), intent(in) :: pid
      real(4), intent(in) :: p0, theta0, phi0
      integer(4), intent(in), optional :: mode

      real(4) p, theta, phi, u, v, w, sum, Kx, Ky, Kz
      integer(4) ix, iy, iz, i, j, k, ilo, ihi, jlo, jhi, klo, khi
      integer(4) xdim, ydim, zdim, mod_
      real(4) plo, pup, dp, thlo, thup, dth, phlo, phup, dph

      if (present(mode)) then
        mod_ = mode
      else
        mod_ = -2   ! use B-spline interpolation
      end if

      getHadesAcceptance = 0.

      if (readHAFTmatrix()==-1) return

      call getDimensions(pid,xdim,ydim,zdim)
      call getLimits(pid,plo,pup,dp,thlo,thup,dth,phlo,phup,dph) 

      if (p0<plo .or. theta0<thlo .or. theta0>thup) return
      if (phi0<phlo) return

      p = min(p0,pup-2.01*dp)  ! level off acceptance at high p
      theta = theta0
      phi = phi0
      if (phi > 60.) phi = mod(phi,60._4)   ! modulo HADES sector

      ix = int(xdim*((p-0.5*dp-plo)/(pup-plo))) + 1      ! floor indices
      iy = int(ydim*((theta-0.5*dth-thlo)/(thup-thlo))) + 1
      iz = int(zdim*((phi-0.5*dph-phlo)/(phup-phlo))) + 1

      select case (mod_)
      case (0,1)  ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      case (2,3,4,-2,-3,-4)
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      case default  ! mode not defined
        return
      end select

      if (ilo<0 .or. jlo<0 .or. klo<0) return
      if (ihi>xdim+1 .or. jhi>ydim+1 .or. khi>zdim+1) return

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



  !*****************************************************************************
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
  !*****************************************************************************
  real(4) function getHadesPairAcceptance(mass0,pt0,rap0,mode)

      real(4), intent(in) :: mass0, pt0, rap0
      integer(4), intent(in), optional :: mode

      real(4) mass, pt, rap, u, v, w, sum, Kx, Ky, Kz
      integer(4) ix, iy, iz, i, j, k, ilo, ihi, jlo, jhi, klo, khi
      integer(4) xdim, ydim, zdim, mod
      real(4) mlo, mup, dm, ptlo, ptup, dpt, raplo, rapup, drap

      getHadesPairAcceptance = 0.

      if (readHAFTPairMatrix()==-1) return

      if (present(mode)) then
        mod = mode  ! (use mode = 0 or 1, otherwise problems at pt=0!)
      else
        mod = 1    ! use trilinear interpolation
      end if

      call getDimensions(51,xdim,ydim,zdim)
      call getLimits(51,mlo,mup,dm,ptlo,ptup,dpt,raplo,rapup,drap) 

      if (mass0<mlo .or. pt0<ptlo .or. pt0>ptup .or. rap0<raplo .or. rap0>rapup) &
        return

      mass = min(mass0,mup-2.01*dm)  ! level off acceptance at high mass
      pt = pt0
      rap = rap0

      ix = int(xdim*((mass-0.5*dm-mlo)/(mup-mlo))) + 1      ! floor indices
      iy = int(ydim*((pt-0.5*dpt-ptlo)/(ptup-ptlo))) + 1
      iz = int(zdim*((rap-0.5*drap-raplo)/(rapup-raplo))) + 1

      select case (mod)
      case (0,1)  ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
        klo = iz
        khi = iz+1
      case (2,3,4,-2,-3,-4)
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
        klo = iz-1
        khi = iz+2
      case default  ! mode not defined
        return
      end select

      if (ilo<0 .or. jlo<0 .or. klo<0) return
      if (ihi>xdim+1 .or. jhi>ydim+1 .or. khi>zdim+1) return

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



  !*****************************************************************************
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
  !*****************************************************************************
  real(4) function kernel(u,mode)

      real(4), intent(in) :: u
      integer(4), intent(in) :: mode

      real(4) ua

      select case (mode)
      case (0)  ! nearest neighbour
        if (u>-0.5 .and. u<=0.5) then
          kernel = 1.
        else 
          kernel = 0.
        end if

      case (1)  ! linear
        ua = abs(u)
        if (ua<1.) then
          kernel = 1.-ua
        else 
          kernel = 0.
        end if

      case (2)  ! quadratic
        ua = abs(u)
        if (ua<=0.5) then
          kernel = 1.0-2.*ua*ua
        else if (ua<=1.5) then
          kernel = 1.5+(ua-2.5)*ua
        else 
          kernel = 0.
        end if

      case (-2)  ! quadratic B-spline
        ua = abs(u)
        if (ua<=0.5) then
          kernel = 0.75-ua*ua
        else if (ua<=1.5) then
          kernel = 1.125+0.5*(ua-3.)*ua
        else 
          kernel = 0.
        end if

      case (3)  ! cubic
        ua = abs(u)
        if (ua<=1.) then
          kernel = 1.+(ua-2.)*ua*ua
        else if (ua<=2.) then
          kernel = 4.+((5.-ua)*ua-8.)*ua
        else 
          kernel = 0.
        end if

      case (-3)  ! cubic B-spline
        ua = abs(u)
        if (ua<=1.) then
          kernel = 2./3.+(0.5*ua-1.)*ua*ua
        else if (ua<=2.) then
          kernel = 4./3.+((1.-ua/6.)*ua-2.)*ua
        else 
          kernel = 0.
        end if

      case (4)  ! cubic Catmull-Rom
        ua = abs(u)
        if (ua<=1.) then
          kernel = 1.+(1.5*ua-2.5)*ua*ua
        else if (ua<=2.) then
          kernel = 2.+((2.5-0.5*ua)*ua-4.)*ua
        else 
          kernel = 0.
        end if

      case (-4)  ! optimal cubic cardinal spline
                 ! (=compromise between blurring and ringing)
        ua = abs(u)
        if (ua<=1.) then
          kernel = 8./9.+(7./6.*ua-2.)*ua*ua
        else if (ua<=2.) then
          kernel = 16./9.+((2.-7./18.*ua)*ua-10./3.)*ua
        else 
          kernel = 0.
        end if

      case default  ! undefined mode
        kernel = 0.

      end select

  end function kernel



  !*****************************************************************************
  !  Opens file in unformatted direct access mode
  !  and reads HADES acceptance matrices (as linearized arrays)
  !*****************************************************************************
  integer(4) function readHAFTmatrix()

      integer(4), parameter :: runit = 77  ! change if input unit is already busy

      integer(4) pid, i, bins
      integer(4) bytes       ! byte counter

      readHAFTmatrix = 0

      if (readflag==1) return

      readflag = 0
      ! set all matrices to 0
      matrix2(:) = 0.
      matrix3(:) = 0.
      matrix8(:) = 0.
      matrix9(:) = 0.
      matrix10(:) = 0.
      matrix12(:) = 0.
      matrix14(:) = 0.

      xdim(:) = 0
      ydim(:) = 0
      zdim(:) = 0
      pmin(:) = 0.
      pmax(:) = 0.
      dp(:) = 0.
      thmin(:) = 0.
      thmax(:) = 0.
      dth(:) = 0.
      phmin(:) = 0.
      phmax(:) = 0.
      dph(:) = 0.

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
      if (pid<0) goto 45
      if (pid<1 .or. pid>nids) goto 50

      ! read acceptance matrix

      bytes = bytes + 4
      read(runit,pos=bytes,err=100) xdim(pid), ydim(pid), zdim(pid)

      bins = xdim(pid)*ydim(pid)*zdim(pid)
      if (bins>size) goto 101  ! check if enough memory available

      bytes = bytes + 3*4
      read(runit,pos=bytes,err=100) pmin(pid),pmax(pid),   &
                                    thmin(pid),thmax(pid), &
                                    phmin(pid),phmax(pid)
      bytes = bytes + 6*4
      if (pid==2) then        ! positron
        read(runit,pos=bytes,err=100) (matrix2(i),i=1,bins)
        write(6,'(''Matrix read for e+'')')
      else if (pid==3) then     ! electron
        read(runit,pos=bytes,err=100) (matrix3(i),i=1,bins)
        write(6,'(''Matrix read for e-'')')
      else if (pid==8) then     ! pi+
        read(runit,pos=bytes,err=100) (matrix8(i),i=1,bins)
        write(6,'(''Matrix read for pi+'')')
      else if (pid==9) then     ! pi-
        read(runit,pos=bytes,err=100) (matrix9(i),i=1,bins)
        write(6,'(''Matrix read for pi-'')')
      else if (pid==10) then    ! K+
        read(runit,pos=bytes,err=100) (matrix10(i),i=1,bins)
        write(6,'(''Matrix read for K+'')')
      else if (pid==12) then    ! K-
        read(runit,pos=bytes,err=100) (matrix12(i),i=1,bins)
        write(6,'(''Matrix read for K-'')')
      else if (pid==14) then    ! proton
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

      ! read resolution parameters

 45   pid = -pid
      if (pid<2 .or. pid>3) goto 102

      bytes = bytes + 4
      read(runit,pos=bytes,err=100) xdimp(pid), ydimp(pid)

      bins = xdimp(pid)*ydimp(pid)
      if (bins>sizep) goto 101  ! check if enough memory available

      bytes = bytes + 2*4
      read(runit,pos=bytes,err=100) pminp(pid),pmaxp(pid), &
                                    thminp(pid),thmaxp(pid)
      bytes = bytes + 4*4
      read(runit,pos=bytes,err=100) ntab ! nb. of parameter tables
      bytes = bytes + 4
      if (pid==2) then          ! positron
        read(runit,pos=bytes,err=100) (par2p1(i),i=1,bins)
        bytes = bytes + bins*4
        if (ntab>1) then
          read(runit,pos=bytes,err=100) (par2p2(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>2) then
          read(runit,pos=bytes,err=100) (par2p3(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>3) then
          read(runit,pos=bytes,err=100) (par2p4(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>4) then
          read(runit,pos=bytes,err=100) (par2p5(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>5) then
          read(runit,pos=bytes,err=100) (par2p6(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        write(6,'(''Resolution tables read for e+'')')
      else if (pid==3) then     ! electron
        read(runit,pos=bytes,err=100) (par3p1(i),i=1,bins)
        bytes = bytes + bins*4
        if (ntab>1) then
          read(runit,pos=bytes,err=100) (par3p2(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>2) then
          read(runit,pos=bytes,err=100) (par3p3(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>3) then
          read(runit,pos=bytes,err=100) (par3p4(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>4) then
          read(runit,pos=bytes,err=100) (par3p5(i),i=1,bins)
          bytes = bytes + bins*4
        end if
        if (ntab>5) then
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

      ! Error opening or reading

 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',trim(fname)
      readHAFTMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',trim(fname)
      readHAFTmatrix = -1
      return
 101  close(runit)
      write(6,*) 'Size error: ', bins, ' >', size,' File = ',trim(fname)
      readHAFTmatrix = -1
      return
 102  close(runit)
      write(6,*) 'PID not yet supported: ', pid, ' File = ',trim(fname)
      readHAFTmatrix = -1
      return
  end function readHAFTmatrix



  !*****************************************************************************
  !  Opens file in unformatted direct access mode
  !  and reads HADES pair acceptance matrix (as linearized array)
  !*****************************************************************************
  integer(4) function readHAFTPairMatrix()

      integer(4), parameter :: runit = 78  ! change if input unit is already busy

      integer(4) i, bins
      integer(4) bytes       ! byte counter

      readHAFTPairMatrix = 0

      if (readflag2==1) return

      readflag2 = 0
      ! set matrix to 0
      matrix51(:) = 0.

      open(unit=runit,file=fname2,access='stream',status='old',err=99)
      bytes=1
      read(runit,pos=bytes,err=100) comment2
      write(6,'(a80)') comment2
      bytes = bytes + 80
      read(runit,pos=bytes,err=100) xdim2, ydim2, zdim2

      bins = xdim2*ydim2*zdim2
      if (bins>size) goto 101  ! check if enough memory available

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

      ! Error opening or reading

 99   close(runit)
      write(6,*) 'Open error on unit ', runit, ' File = ',trim(fname2)
      readHAFTPairMatrix = -1
      return
 100  close(runit)
      write(6,*) 'Read error on unit ', runit, ' File = ',trim(fname2)
      readHAFTPairMatrix = -1
      return
 101  close(runit)
      write(6,*) 'Size error: ', bins,' >',size, ' File = ',trim(fname2)
      readHAFTPairMatrix = -1
      return
  end function readHAFTPairMatrix



  !*****************************************************************************
  !  Sets name of input file containing the filter
  !*****************************************************************************
  subroutine setFileName(name)
      character*(*), intent(in) :: name

      integer(4) dummy

      fname = name
      dummy = readHAFTmatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname |'',a80,''|'')') fname
  end subroutine setFileName



  !*****************************************************************************
  !  Sets name of input file containing the filter
  !*****************************************************************************
  subroutine setPairFileName(name)
      character*(*), intent(in) :: name

      integer(4) dummy

      fname2 = name
      dummy = readHAFTPairMatrix()

!      write(6,'(''name  |'',a80,''|'')') name
!      write(6,'(''fname2 |'',a80,''|'')') fname2
  end subroutine setPairFileName



  !*****************************************************************************
  !  Returns acceptance value at cell (i,j,k) of linearized matrix
  !  for particle ID
  !*****************************************************************************
  real(4) function getMatrixVal(i,j,k,pid)
      integer(4), intent(in) :: i, j, k, pid

      integer(4) xdi, ydi, zdi, i1, j1, k1, ilin

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



  !*****************************************************************************
  !  Return the lower and upper limits, and step sizes of the table
  !*****************************************************************************
  subroutine getLimits(pid,xlo,xhi,dx,ylo,yhi,dy,zlo,zhi,dz) 
      integer(4), intent(in) :: pid
      real(4), intent(out) :: xlo, xhi, dx, ylo, yhi, dy, zlo, zhi, dz

      if (pid==51) then  ! pair
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



  !*****************************************************************************
  !  Return the dimensions of a table of particle pid
  !*****************************************************************************
  subroutine getDimensions(pid,nx,ny,nz)
      integer(4), intent(in)  :: pid
      integer(4), intent(out) :: nx, ny, nz

      if (pid==51) then  ! pair
        nx = xdim2
        ny = ydim2
        nz = zdim2
      else
        nx = xdim(pid)
        ny = ydim(pid)
        nz = zdim(pid)
      end if
  end subroutine getDimensions



  !*****************************************************************************
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
  !*****************************************************************************
  subroutine smearHades4Momentum(mom,mode,pid)

      real(4), intent(inout) :: mom(4)
      integer(4), intent(in) :: mode,pid

      integer(4) i
      real(4) mass, mass2, pt, pt2, ptot, ptot2, theta, phi, sinth
      real(4) sigp, betainv, sigms, sigms2, sigthms, sigphms, ploss, respar(10)
      real(4), parameter :: r2d = 180./pi

      if (readflag==0) then
        if (readHAFTmatrix()==-1) return
      end if

      pt2 = mom(1)**2 + mom(2)**2
      pt = sqrt(pt2)
      ptot2 = pt2 + mom(3)**2
      ptot = sqrt(ptot2)                      ! total momentum in GeV/c
      mass2 = mom(4)**2 - ptot2               ! particle mass does not change
      if (mass2<0.) mass2 = 0.
      mass = sqrt(mass2)
      betainv = 1.
      if (mass2>0.) betainv = sqrt(1. + mass2/ptot2)  ! 1/beta
      sigms = 0.0136*betainv/ptot*sqrt(XX0)*(1+0.038*log(XX0))
                                              ! multiple scattering angle
      sigms2 = sigms*sigms
      sigthms = sqrt(sigth*sigth+sigms2)      ! add quadratically to resolution
      sigphms = sqrt(sigph*sigph+sigms2)

      if (pid==2 .or. pid==3) then  ! leptons, use Ar+KCl embedding values
         sigthms = 0.00246/ptot + 0.00045
         sigphms = 0.00294/ptot + 0.00059
      endif

      theta = acos(mom(3)/ptot)               ! polar angle in radian
      sinth = sin(theta)
      phi = 0.
      if (pt>0.) phi = acos(mom(1)/pt)
      if (mom(2)<0.) phi = twopi - phi     ! azimuthal angle in radian

      !  If resolution parameters are available, use dedicated smearing

      if (resflag(pid)==1) then   ! resolution parameters are loaded for pid?

!        write(6,*) mom(1), mom(2), mom(3), mass, ptot, r2d*theta
         do i=1,ntab    ! look up parameters from tables
            respar(i) = param(ptot,r2d*theta,pid,i)
!           write(6,*) i, respar(i)
         end do

         ptot = ptot*(1.+sampleMP(respar,2._4))  ! randomize momentum

      else  ! use default gaussian smearing

         if (mode<1 .or. mode>3) return  ! unknown mode
         sigp = 0.01*ptot*(sigpA(mode)+sigpB(mode)*ptot)  ! momentum resolution
         ptot = max(0.,sampleGauss(ptot,sigp))  ! smear total momentum

         if (pid==2 .or. pid==3) then
            ploss = 0.0018
         else
            ploss = 0.0
         end if
         if (ptot>0.1) then
            if (pid==2) then
               ploss = ploss + 0.005*(1.-exp(-(ptot-0.1)/0.053)) ! positron
            else if (pid==3) then
               ploss = ploss + 0.006*(1.-exp(-(ptot-0.1)/0.078)) ! electron
            end if
         end if
         ptot = ptot - ploss                    ! Eloss of electrons in target

      end if

      theta = abs(sampleGauss(theta,sigthms)) ! smear polar angle
      if (theta>pi) theta = twopi - theta  ! if > pi, mirror angle
      if (sinth>0.01) phi = sampleGauss(phi,sigphms/sinth) ! smear azimuth
      if (phi<0.) phi = phi + twopi        ! and check if within range 
      if (phi>twopi) phi = phi - twopi

      sinth = sin(theta)
      mom(1) = ptot*sinth*cos(phi)            ! new momentum components
      mom(2) = ptot*sinth*sin(phi)
      mom(3) = ptot*cos(theta)
      mom(4) = sqrt(sum(mom(1:3)**2) + mass2)  ! total energy

  end subroutine smearHades4Momentum



  subroutine smearHadesMomentum(p,mode,pid)
      real, intent(inout) :: p(0:3)
      integer(4), intent(in) :: mode, pid

      real(4) mom4(4)

      mom4(1:3) = p(1:3)
      mom4(4) = p(0)

      call smearHades4Momentum(mom4,mode,pid)

      p(1:3) = mom4(1:3)
      p(0) = mom4(4)
  end subroutine smearHadesMomentum



  !*****************************************************************************
  !  Apply Hades momentum resolution to a 3-momentum (calculate multiple
  !  scattering assuming the particle is an electron)
  !*****************************************************************************
  subroutine smearHades3Momentum(mom3,mode,pid)
      real(4), intent(inout) :: mom3(3)
      integer(4), intent(in) :: mode,pid

      real(4) mom4(4), mass

      if (mode<1 .or. mode>3) return    ! unknown mode

      if (pid==2 .or. pid==3) then        ! e+ or e-
         mass = 0.000510998918
      else if (pid==8 .or. pid==9) then   ! pi+ or pi- 
         mass = 0.13957018
      else if (pid==11 .or. pid==12) then ! K+ or K- 
         mass = 0.493677
      else if (pid==14) then                ! proton
         mass = 0.938272029
      else
         return  ! particle is not supported
      end if

      mom4(1:3) = mom3(1:3)
      mom4(4) = sqrt(sum(mom3(:)**2) + mass**2)
      call smearHades4Momentum(mom4,mode,pid)
      mom3(1:3) = mom4(1:3)
  end subroutine smearHades3Momentum



  !*****************************************************************************
  !  Apply Hades momentum resolution to a pair (calculate multiple
  !  scattering assuming the particle is an electron)
  !*****************************************************************************
  subroutine smearHadesPair(pair,mode)
      real(4), intent(inout) :: pair(3)
      integer(4), intent(in) :: mode

      real(4) m, pt, rap, sigpt, sigm, sigrap, respar(10)

      real(4), parameter :: mtab(10) = (/ 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2 /)  ! mass grid
      real(4), parameter :: par1(10) = (/ 0.0077, -0.0082, -0.0125, -0.0120, -0.0114, &
                                        -0.0106, -0.0098, -0.0085, -0.0078, -0.0075 /)
      real(4), parameter :: par2(10) = (/ 0.0820, 0.0460, 0.0260, 0.0210, 0.0190, &
                                         0.0183, 0.0182, 0.0181, 0.0180, 0.0180 /)
      real(4), parameter :: par3(10) = (/ 13.5, 16.2, 19.9, 20.2, 19.2, &
                                         18.0, 16.9, 14.8, 12.8, 11.0 /)
      real(4), parameter :: par4(10) = (/-10.0, -15.4, -20.4, -23.1, -22.3, &
                                        -21.6, -20.6, -19.4, -18.2, -17.9 /)
      real(4), parameter :: par5(10) = (/ 18.1, 11.6, 10.4, 10.0, 9.4, 8.5, 7.8, 7.0, 6.2, 5.7 /)

      if (readflag==0) then
        if (readHAFTmatrix()==-1) return
      end if

      m = pair(1)
      pt = pair(2)
      rap = pair(3)

      if (mode==4) then  ! use skewd mass function
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


  !*****************************************************************************
  !  Return random number according to a normal distribution.
  !*****************************************************************************
  real(4) function sampleGauss(mean,sigma)
      real(4), intent(in) :: mean, sigma

      real(4) theta, r(2)

      sampleGauss = mean
      if (sigma<=0.) return
      call random_number(r)
      theta = twopi*r(1)
      sampleGauss = mean + sigma*cos(theta)*sqrt(-2.*log(r(2)))
  end function sampleGauss



  !*****************************************************************************
  !     HADES momentum spread (pRec-pSim)/pSim
  !*****************************************************************************
  real(4) function momSpread(x,respar,ns)

      real(4), intent(in) :: x, respar(10), ns

      real(4) pos, sig, left, right, farleft, argn, argp, argn2, e2, amp

      e2 = exp(-0.5*ns*ns)

      pos = respar(1)     ! Mean
      sig = respar(2)     ! Sigma
      left = respar(3)    ! Par3 (>0)
      right = respar(4)   ! Par4 (<0)
      farleft = respar(5) ! Par5 (>0)

      if (x>=(pos-ns*sig) .or. x<=(-lg10/left+pos-ns*sig)) then
         argn = 0.
      else 
         argn = 1.
      end if

      if (x>=(pos-ns*sig)) then
        argp = 1.
      else
        argp = 0.
      end if

      if (x>(-lg10/left+pos-ns*sig)) then
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



  !*****************************************************************************
  !  Return random number according to the normalized HADES momentum distribution.
  !*****************************************************************************
  real(4) function sampleMP(respar,ns)

      real(4), intent(in) :: respar(10), ns

      real(4) pos, sig, left, right, farleft, A(0:3), F(0:3)
      real(4) dx, ftest, r1, r2, r3, e2
      integer(4) cnt, cnt1, cnt2, cnt3

      e2 = exp(-0.5*ns*ns)

      pos = respar(1)      ! centroid
      sig = respar(2)      ! width
      left = respar(3)     ! left slope
      right = respar(4)    ! right slope
      farleft = respar(5)  ! far left slope

      ! compute function amplitudes
      A(0) = 0.1*(1. + e2)
      A(1) = 1. + e2
      A(2) = 1. + exp(right*ns*sig)
      A(3) = (e2 + exp(right*2.*ns*sig))/exp(right*2.*ns*sig)

      ! compute function areas
      F(0) = A(0)/farleft * (1. - exp(farleft*(lg10/left+ns*sig-pos-1.)))  ! [-1, 1/10left]
      F(1) = A(1)/left * 9./10.                                            ! [1/10, pos-ns*sig]
      F(2) = A(2) * 2.*ns*sig                                              ! [pos-ns*sig, pos+ns*sig]
      F(3) = A(3)/right * (exp(right*(1.-pos+ns*sig))-exp(right*2.*ns*sig))! [pos + ns*sig, 1]

      F(0:3) = F(0:3)/sum(F)   ! normalize areas

      ! sample dx by comparing with piece-wise function

      do cnt=1,1000    ! allow max 1000 trials
         cnt1 = 0
         cnt2 = 0
         cnt3 = 0
         call random_number(r1)
         ! select region and sample test function
         if  (r1 < F(0)) then             ! far left tail

            do
               cnt1 = cnt1 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/farleft - ns*sig + pos - lg10/left
               if (cnt1==1000) write(6,*) 'cnt1=1000 ', pos, sig, farleft
               if (dx>=-1. .or. cnt1>=1000) exit  ! limit range to >=-1
            end do
            ftest = A(0) * exp(farleft*(dx+lg10/left-pos+ns*sig))

         else if (r1 < sum(F(0:1))) then      ! left tail

            do
               cnt2 = cnt2 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/left - ns*sig + pos
               if (cnt2==1000) write(6,*) 'cnt2=1000', pos, sig, left
               if (dx>=-lg10/left-ns*sig+pos .or. cnt2>=1000) exit
            end do
            ftest = A(1) * exp(left*(dx-pos+ns*sig))

         else if (r1 < sum(F(0:2))) then   ! peak region

            call random_number(r2)
            r2 = r2 - 0.5
            dx = 2.*ns*sig*r2 + pos
            ftest = A(2)

         else                            ! right tail

            do
               cnt3 = cnt3 + 1
               call random_number(r2)
               r2 = log(r2)
               dx = r2/right + ns*sig + pos
               if (cnt3==1000) write(6,*) 'cnt3=1000', pos, sig, right
               if (dx<=1. .or. cnt3>=1000) exit  ! limit range to <=1
            end do
            ftest = A(3) * exp(right*(dx-pos+ns*sig))

         end if

         ! do rejection test
         sampleMP = dx

         call random_number(r3)
         if ( r3<momSpread(dx,respar,ns)/ftest ) return

      end do

      write(6,*) 'cnt=1000'
      sampleMP = 0.

  end function sampleMP



  !*****************************************************************************
  !     Set momentum resolution parameters
  !*****************************************************************************
  subroutine setResolutionParameters(mode,a,b)
      real(4), intent(in) :: a, b
      integer(4), intent(in) :: mode

      if (mode>=1 .and. mode<=3) then
         sigpA(mode) = a
         sigpB(mode) = b
      endif

      write(6,*) 'mode=1   ', ' ', sigpA(1), ' ', sigpB(1)
      write(6,*) 'mode=2   ', ' ', sigpA(2), ' ', sigpB(2)
      write(6,*) 'mode=3   ', ' ', sigpA(3), ' ', sigpB(3)
  end subroutine setResolutionParameters



  !*****************************************************************************
  !     Set angular resolution parameters
  !*****************************************************************************
  subroutine setAngularResolutionParameters(stheta,sphi)
      real(4), intent(in) :: stheta, sphi

      if (stheta>0. .and. sphi>0.) then
         sigth = stheta
         sigph = sphi
      endif

      write(6,*) 'dTh, dPh:', ' ', sigth, ' ', sigph
  end subroutine setAngularResolutionParameters



  !*****************************************************************************
  !     Interpolate resolution parameter table as function
  !     of momentum and theta (pin in GeV/c and theta in degree)
  !*****************************************************************************
  real(4) function param(pin,thin,pid,itab)
      real(4), intent(in) :: pin, thin
      integer(4), intent(in) :: pid, itab

      integer(4) xdi, ydi, i, j, ix, iy, ilo, ihi, jlo, jhi, mod
      real(4) p, th, plo, pup, dp0, thlo, thup, dth0
      real(4) sum, u, v, Kx, Ky

      mod = 1
      param = 0.
      if (pid<1 .or. pid>nids) return

      p = pin
      th = thin
      plo = pminp(pid)
      pup = pmaxp(pid)
      dp0 = dpp(pid)
      thlo = thminp(pid)
      thup = thmaxp(pid)
      dth0 = dthp(pid)
      if (p<plo) p = plo   ! safety fence
      if (p>pup) p = pup
      if (th<thlo) th = thlo
      if (th>thup) th = thup

      xdi = xdimp(pid) ! get table dimensions
      ydi = ydimp(pid)

      ix = int(xdi*((p-0.5*dp0-plo)/(pup-plo))) + 1      ! floor indices
      iy = int(ydi*((th-0.5*dth0-thlo)/(thup-thlo))) + 1

      if (mod==0 .or. mod==1) then   ! set summation limits
        ilo = ix
        ihi = ix+1
        jlo = iy
        jhi = iy+1
      else if (abs(mod)==2 .or. abs(mod)==3 .or. abs(mod)==4) then
        ilo = ix-1
        ihi = ix+2
        jlo = iy-1
        jhi = iy+2
      else  ! mode not defined
        return
      end if

      if (ilo<0 .or. jlo<0) return
      if (ihi>xdi+1 .or. jhi>ydi+1) return

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



  !*****************************************************************************
  !  Returns acceptance value at cell (i,j) of linearized
  !  parameter table for particle ID
  !*****************************************************************************
  real(4) function getTableVal(i,j,pid,itab)
      integer(4), intent(in) :: i, j, pid, itab

      integer(4) xdi, ydi, i1, j1, ilin

      getTableVal = 0.
      if (pid<1 .or. pid>nids) return

      xdi = xdimp(pid) ! get table dimensions
      ydi = ydimp(pid)

      i1 = min(max(1,i),xdi)  ! Make sure indexes stay within table.
      j1 = min(max(1,j),ydi)  ! This effectively extrapolates matrix

      ilin = i1+xdi*(j1-1)    ! linearized index

      if (pid==2) then          ! positron
        if (itab==1) getTableVal = par2p1(ilin)
        if (itab==2) getTableVal = par2p2(ilin)
        if (itab==3) getTableVal = par2p3(ilin)
        if (itab==4) getTableVal = par2p4(ilin)
        if (itab==5) getTableVal = par2p5(ilin)
        if (itab==6) getTableVal = par2p6(ilin)
      else if (pid==3) then     ! electron
        if (itab==1) getTableVal = par3p1(ilin)
        if (itab==2) getTableVal = par3p2(ilin)
        if (itab==3) getTableVal = par3p3(ilin)
        if (itab==4) getTableVal = par3p4(ilin)
        if (itab==5) getTableVal = par3p5(ilin)
        if (itab==6) getTableVal = par3p6(ilin)
      end if
  end function getTableVal



  !*****************************************************************************
  !     linear interpolation in table (xtab,ytab)
  !*****************************************************************************
  real(4) function interpol(x,xtab,ytab,n)
      real(4), intent(in) :: x, xtab(*), ytab(*)
      integer(4), intent(in) :: n

      integer(4) i
      real(4) a, b

      if (x<=xtab(1)) then ! below table range
        interpol = ytab(1)
        return
      else if (x>=xtab(n)) then ! above table range
        interpol = ytab(n)
        return
      end if

      do i=2,n
        interpol = ytab(i)
        if (x==xtab(i)) return
        if (x<xtab(i)) exit
      end do

      a = ytab(i-1)
      b = (ytab(i)-ytab(i-1))/(xtab(i)-xtab(i-1))
      interpol = a + (x-xtab(i-1))*b
  end function interpol


end module
