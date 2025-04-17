      program SOLVER33_TEST_SMP

      use solver33
      use hpcmw_all
      implicit REAL*8(A-H,O-Z)
      integer, dimension(:), allocatable :: OLDtoNEWpe
!C
!C +-------+
!C | INIT. |
!C +-------+
!C=== 
      call HPCMW_INIT
      call INPUT_CNTL

      allocate (OLDtoNEWpe(PETOT))

      do ITERkk= 1, 1
        call MPI_BARRIER (MPI_COMM_WORLD, errno)
        if (ITERkk.eq. 1) NCOLORk= 0

      call INPUT_GRID (OLDtoNEWpe, ITERkk)
!C===

!C
!C +---------------------+
!C | matrix connectivity |
!C +---------------------+
!C===
      call MAT_CON0 (ITERkk)
!C===

!C
!C +-------------+
!C | RE-ORDERING |
!C +-------------+
!C===
      call MAT_CON1 
!C===
!C
!C +-----------------+
!C | MATRIX assemble |
!C +-----------------+h
!C===
      call MAT_ASS_MAIN
!C===

!C
!C +---------------------+
!C | BOUNDARY CONDITIONs |
!C +---------------------+
!C===
      call MAT_ASS_BC 
!C===

!C
!C +--------+
!C | SOLVER |
!C +--------+
!C===
      call SOLVE33 (hpcmwIarray, hpcmwRarray, ITERkk)

      if (my_rank.eq.PETOT-1) then
        i= N
        write (*,'(i8,3(1pe16.6))') i,X(3*i-2),X(3*i-1),X(3*i)
      endif

      enddo

      call HPCMW_FINALIZE
      end program SOLVER33_TEST_SMP

!C
!C***
!C*** ERROR_EXIT
!C***
!C
      subroutine ERROR_EXIT (IFLAG, my_rank)
      return
      end
      

