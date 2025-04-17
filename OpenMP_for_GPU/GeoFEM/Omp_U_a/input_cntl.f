!======================================================================!
!                                                                      !
!  Software Name : HPCMW-SOLVER-SMP-TEST Ver.1.0                       !
!                                                                      !
!    main program:                                                     !
!    subroutine  : INPUT_CNTL                                          !
!    module name :                                                     !
!    class name  :                                                     !
!                                                                      !
!    coded by Kengo Nakajima   (RIST)   2003/03/27                     !
!                                                                      !
!    Contact address :  The University of Tokyo, FSIS project          !
!                                                                      !
!    "High-Performance Computing Middleware (HPC-MW)" Group.           !
!                                                                      !
!======================================================================!
!C
!C***
!C*** INPUT_CNTL
!C***
!C
      subroutine INPUT_CNTL
      use hpcmw_solver_cntl
      use hpcmw_solver_matrix
      use hpcmw_fem_cntl

      implicit REAL*8 (A-H,O-Z)

      PEsmpTOT= 16
      FTflag  = 1

      METHOD  =    1
      PRECOND =    0

!      open (11,file='COLOR.DAT',status='unknown')
!        read (11,*) NCOLORk
!      close (11)

      NCOLORk= -0

      iterPREmax=  2
      ITER      = 2000000

      SIGMA_DIAG= 1.d0
      SIGMA     = 0.d0
      RESID     = 1.d-8
      NSET      = 0

      hpcmwRarray(1)= RESID
      hpcmwRarray(2)= SIGMA_DIAG
      hpcmwRarray(3)= SIGMA

      hpcmwIarray(1)= ITER
      hpcmwIarray(2)= METHOD
      hpcmwIarray(3)= PRECOND
      hpcmwIarray(4)= NSET
      hpcmwIarray(5)= iterPREmax

      return
      end
