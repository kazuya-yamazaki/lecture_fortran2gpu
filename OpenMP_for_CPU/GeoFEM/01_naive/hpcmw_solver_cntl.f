!======================================================================!
!                                                                      !
!  Software Name : HPCMW-SOLVER-SMP-TEST Ver.1.0                       !
!                                                                      !
!    main program:                                                     !
!    subroutine  :                                                     !
!    module name : hpcmw_solver_cntl                                   !
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
!C*** HPCmw_solver_cntl
!C***
!C
      module hpcmw_solver_cntl
      use hpcmw_util

      integer(kind=kint) :: METHOD, PRECOND, NSET, iterPREmax
      integer(kind=kint) :: ITER, ITERactual

      real   (kind=kreal) :: RESID, SIGMA_DIAG, SIGMA

      end module hpcmw_solver_cntl
