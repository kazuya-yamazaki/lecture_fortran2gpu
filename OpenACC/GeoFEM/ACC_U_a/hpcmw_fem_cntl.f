!======================================================================!
!                                                                      !
!  Software Name : HPCMW-SOLVER-SMP-TEST Ver.1.0                       !
!                                                                      !
!    main program:                                                     !
!    module name : hpcmw_fem_cntl                                      !
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
!C*** HPCmw_fem_cntl
!C***
!C
      module hpcmw_fem_cntl
      use hpcmw_util

      integer(kind=kint ), dimension(100) :: hpcmwIarray
      real   (kind=kreal), dimension(100) :: hpcmwRarray

      real   (kind=kreal), parameter :: O8th= 0.125d0

      end module hpcmw_fem_cntl
