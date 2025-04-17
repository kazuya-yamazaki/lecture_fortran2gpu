!======================================================================!
!                                                                      !
!  Software Name : HPCMW-SOLVER-SMP-TEST Ver.1.0                       !
!                                                                      !
!    main program:                                                     !
!    subroutine  :                                                     !
!    module name : hpcmw_util                                          !
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
!C*** HPCmw_util
!C***
!C
!C    MPI settings
!C    FILE names
!C
      module hpcmw_util
      include 'mpif.h'
      include 'precision.inc'

      integer :: PEsmpTOT
      integer :: PETOT, my_rank, SOLVER_COMM, errno

      character(len= 4 ) ::  penum, penum_left
      character(len= 9 ) ::  fname, rname

      end module hpcmw_util
