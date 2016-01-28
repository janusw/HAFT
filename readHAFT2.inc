c  HAFT declaration of acceptance matrix arrays and resolution tables
c
c  Code version 2.0 of February 15, 2011
c
c
c  The dimensions MUST match all array sizes in the file!
c
      integer*4 size, sizep, nids
      parameter (size=250000) ! <<== change if < xdim*ydim*zdim
      parameter (sizep=1000)  ! <<== change if < xdimp*ydimp
      parameter (nids=14)     ! <<== change if < max id
c
      character*80 fname, fname2
      character*80 comment, comment2
      integer*4 xdim(nids), ydim(nids), zdim(nids)
      real*4 dp(nids), dth(nids), dph(nids)
      real*4 pmin(nids), pmax(nids), thmin(nids),
     +       thmax(nids), phmin(nids), phmax(nids)
      real*4 matrix2(size) ! matrices are declared for
      real*4 matrix3(size) ! e+, e-, pi+, pi-, K+, K- and p
      real*4 matrix8(size)
      real*4 matrix9(size)
      real*4 matrix10(size)
      real*4 matrix12(size)
      real*4 matrix14(size)
      integer*4 iseed
      integer*4 readflag, readflag2
      real*4 sigpA(3), sigpB(3), sigth, sigph, XX0
      integer*4 xdim2, ydim2, zdim2
      real*4 dm, dpt, drap
      real*4 mmin, mmax, ptmin, ptmax, rapmin, rapmax
      real*4 matrix51(size)

      integer*4 resflag(nids)
      integer*4 xdimp(nids), ydimp(nids)
      real*4 dpp(nids), dthp(nids)
      real*4 pminp(nids), pmaxp(nids), thminp(nids),
     +       thmaxp(nids)
      integer*4 ntab
      real*4 par2p1(sizep), par3p1(sizep) 
      real*4 par2p2(sizep), par3p2(sizep)
      real*4 par2p3(sizep), par3p3(sizep)
      real*4 par2p4(sizep), par3p4(sizep)
      real*4 par2p5(sizep), par3p5(sizep)
      real*4 par2p6(sizep), par3p6(sizep)
      common /haft/ fname, comment, xdim, ydim, zdim, dp, dth, dph,
     +              pmin, pmax, thmin, thmax, phmin, phmax,
     +              matrix2, matrix3, matrix8, matrix9,
     +              matrix10, matrix12, matrix14,
     +              iseed, sigpA, sigpB, sigth, sigph, XX0,
     +              fname2, comment2, xdim2, ydim2, zdim2,
     +              dm, dpt, drap,
     +              mmin, mmax, ptmin, ptmax, rapmin, rapmax,
     +              matrix51, readflag, readflag2,
     +              ntab, resflag, xdimp, ydimp, dpp, dthp,
     +              pminp, pmaxp, thminp, thmaxp,
     +              par2p1, par2p2, par2p3, par2p4, par2p5, par2p6,
     +              par3p1, par3p2, par3p3, par3p4, par3p5, par3p6