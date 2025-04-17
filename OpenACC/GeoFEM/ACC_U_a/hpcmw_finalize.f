!======================================================================!
!                                                                      !
!  Software Name : HPCMW-SOLVER-SMP-TEST Ver.1.0                       !
!                                                                      !
!    main program:                                                     !
!    subroutine  : hpcmw_finalize                                      !
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
!C*** HPCMW_FINALIZE
!C***
!C
      subroutine HPCMW_FINALIZE
      use hpcmw_fem_cntl
      implicit REAL*8 (A-H,O-Z)

      call MPI_FINALIZE (ierr)
      if (my_rank.eq.0) stop ' * normal termination'

      return
      end
