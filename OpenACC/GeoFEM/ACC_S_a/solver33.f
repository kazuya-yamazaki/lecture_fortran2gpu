      module SOLVER33

      contains
        subroutine SOLVE33 (hpcmwIarray, hpcmwRarray, ITERkk) 

        use hpcmw_solver_matrix
        use hpcmw_solver_cntl
        use hpcmw_fem_mesh

        use solver_CG_3_SMP_novec

        implicit REAL*8 (A-H,O-Z)

        real(kind=kreal), dimension(3,3) :: ALU
        real(kind=kreal), dimension(3)   :: PW  

        integer :: ERROR, ICFLAG
        character(len=char_length) :: BUF

        data ICFLAG/0/

        integer(kind=kint ), dimension(:) :: hpcmwIarray
        real   (kind=kreal), dimension(:) :: hpcmwRarray


!C
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER      = hpcmwIarray(1)
        METHOD    = hpcmwIarray(2)
        PRECOND   = hpcmwIarray(3)
        NSET      = hpcmwIarray(4)
        iterPREmax= hpcmwIarray(5)

        RESID     = hpcmwRarray(1)
        SIGMA_DIAG= hpcmwRarray(2)

        if (iterPREmax.lt.1) iterPREmax= 1
        if (iterPREmax.gt.4) iterPREmax= 4
!C===

!C
!C +-----------+
!C | BLOCK LUs |
!C +-----------+
!C===
!$acc kernels present(D,ALUG)
!$acc loop independent private(ALU,PW) 
          do ii= 1, N
            ALU(1,1)= D(9*ii-8)*SIGMA_DIAG
            ALU(1,2)= D(9*ii-7)
            ALU(1,3)= D(9*ii-6)
            ALU(2,1)= D(9*ii-5)
            ALU(2,2)= D(9*ii-4)*SIGMA_DIAG
            ALU(2,3)= D(9*ii-3)
            ALU(3,1)= D(9*ii-2)
            ALU(3,2)= D(9*ii-1)
            ALU(3,3)= D(9*ii  )*SIGMA_DIAG
            do k= 1, 3
               L = k
              ALO= dabs(ALU(L,k))
              do i= k+1, 3
                if (dabs(ALU(i,k)).gt.ALO) then
                   L = i
                  ALO= dabs(ALU(L,k))
                endif
              enddo

              ALU(k,k)= 1.d0/ALU(k,k)
              do i= k+1, 3
                ALU(i,k)= ALU(i,k) * ALU(k,k)
                do j= k+1, 3
                  PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
                enddo
                do j= k+1, 3
                  ALU(i,j)= PW(j)
                enddo
              enddo
            enddo
            ALUG(9*ii-8)= ALU(1,1)
            ALUG(9*ii-7)= ALU(1,2)
            ALUG(9*ii-6)= ALU(1,3)
            ALUG(9*ii-5)= ALU(2,1)
            ALUG(9*ii-4)= ALU(2,2)
            ALUG(9*ii-3)= ALU(2,3)
            ALUG(9*ii-2)= ALU(3,1)
            ALUG(9*ii-1)= ALU(3,2)
            ALUG(9*ii  )= ALU(3,3)
          enddo
!$acc end kernels
!C===

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
        call CG_3_SMP_novec                                             &
     &   ( N, NP, NPL, NPU, PEsmpTOT, NHYP, IVECT,                      &
     &     STACKmc,  NEWtoOLD, OLDtoNEW,                                &
     &     D, AL, indexL, itemL, AU, indexU, itemU,                     &
     &     B, X,  ALUG, WW, WS, WR, RESID, ITER, ERROR,                 &
     &           my_rank, NEIBPETOT, NEIBPE,                            &
     &           NOD_STACK_IMPORT, NOD_IMPORT,                          &
     &           NOD_STACK_EXPORT, NOD_EXPORT,                          &
     &           SOLVER_COMM , PRECOND, iterPREmax)     
      ITERactual= ITER
!C===

      end subroutine SOLVE33
      end module SOLVER33
