!C
!C*** 
!C*** module solver_CG_3_SMP_novec
!C***
!C
      module solver_CG_3_SMP_novec
      contains
!C
!C*** CG_3_SMP_novec
!C
!C    CG_3_SMP_novec solves the linear system Ax = b with 3*3 block matrix 
!C    using the Conjugate Gradient iterative method with the following
!C    preconditioners for SMP nodes:
!C
      subroutine CG_3_SMP_novec                                         &
     &                 (N, NP, NPL, NPU, PEsmpTOT, NVECT, IVECT,        &
     &                  STACKmc, NtoO, OtoN,                            &
     &                  D, AL, INL, IAL, AU, INU, IAU,                  &
     &                  B,  X, ALU, WW, WS, WR, RESID,  ITER, ERROR,    &
     &                  my_rank, NEIBPETOT, NEIBPE,                     &
     &                  STACK_IMPORT, NOD_IMPORT,                       &
     &                  STACK_EXPORT, NOD_EXPORT,                       &
     &                  SOLVER_COMM , PRECOND, iterPREmax)     

      use  solver_SR_3

      implicit REAL*8(A-H,O-Z)
      include  'precision.inc'
      include  'mpif.h'

      integer(kind=kint ), intent(in):: N, NP, NPU, NPL, NVECT, my_rank
      integer(kind=kint ), intent(in):: PEsmpTOT
      integer(kind=kint ), intent(in):: NEIBPETOT
      integer(kind=kint ), intent(in):: SOLVER_COMM
      integer(kind=kint ), intent(in):: PRECOND

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID

      real(kind=kreal), dimension(3*NP )   , intent(inout):: B, X
      real(kind=kreal), dimension(9*NPL), intent(inout):: AL
      real(kind=kreal), dimension(9*NPU), intent(inout):: AU
      real(kind=kreal), dimension(9*NP), intent(inout):: ALU
      real(kind=kreal), dimension(9*NP), intent(inout):: D

      integer(kind=kint ), dimension(0:NP ) ,intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU
      integer(kind=kint ), dimension(0:PEsmpTOT,NVECT),intent(in) 
     &                                                  :: STACKmc  

      real(kind=kreal), dimension(3*NP,4), intent(inout):: WW
      real(kind=kreal), dimension(3*NP)  , intent(inout):: WS, WR

      integer(kind=kint ):: NEIBPE(26)
      integer(kind=kint ):: STACK_IMPORT(0:26), NOD_IMPORT(NP)
      integer(kind=kint ):: STACK_EXPORT(0:26), NOD_EXPORT(NP)

      integer(kind=kint), dimension(NP)   :: NtoO, OtoN
      integer(kind=kint), dimension(0:NP) :: IVECT

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: ZP= 4

      integer(kind=kint ) :: MAXIT, IFLAG, iterPREmax

      real   (kind=kreal) :: TOL, W, SS
      data IFLAG/0/

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0

      COMMtime= 0.d0
      COMPtime= 0.d0

      MAXIT  = ITER
       TOL   = RESID           

      istart = STACK_EXPORT(0) + 1
      iend = STACK_EXPORT(NEIBPETOT)

!$omp parallel do private(in)
      do i= 1, NP
        in= OtoN(i)
        WS(3*in-2)= B(3*i-2) 
        WS(3*in-1)= B(3*i-1) 
        WS(3*in  )= B(3*i  ) 
      enddo
!$omp end parallel do 

!$omp parallel do 
      do i= 1, NP
        B(3*i-2)= WS(3*i-2) 
        B(3*i-1)= WS(3*i-1) 
        B(3*i  )= WS(3*i  ) 
      enddo
!$omp end parallel do 

!$omp parallel do 
      do i= 1, NP
        WW(3*i-2,1)= 0.d0
        WW(3*i-2,2)= 0.d0
        WW(3*i-2,3)= 0.d0
        WW(3*i-2,4)= 0.d0
        WW(3*i-1,1)= 0.d0
        WW(3*i-1,2)= 0.d0
        WW(3*i-1,3)= 0.d0
        WW(3*i-1,4)= 0.d0
        WW(3*i  ,1)= 0.d0
        WW(3*i  ,2)= 0.d0
        WW(3*i  ,3)= 0.d0
        WW(3*i  ,4)= 0.d0
        WS(3*i-2)= 0.d0
        WS(3*i-1)= 0.d0
        WS(3*i  )= 0.d0
        WR(3*i-2)= 0.d0
        WR(3*i-1)= 0.d0
        WR(3*i  )= 0.d0
        X (3*i-2)= 0.d0
        X (3*i-1)= 0.d0
        X (3*i  )= 0.d0
      enddo
!$omp end parallel do 
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

      call SOLVER_SEND_RECV_3                                           &
     &   ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,          &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)
!C
!C-- BEGIN calculation

      do jv= 1, NVECT
        jS= IVECT(jv-1) + 1
        jE= IVECT(jv)

      do ip= 1,PEsmpTOT
        jS= STACKmc(ip-1,jv) + 1
        jE= STACKmc(ip  ,jv)
!$omp parallel do private(X1,X2,X3,WVAL1,WVAL2,WVAL3,k,i)
      do j= jS, jE
           X1= X(3*j-2)
           X2= X(3*j-1)
           X3= X(3*j  )
        WVAL1= B(3*j-2) - D(9*j-8)*X1 - D(9*j-7)*X2 - D(9*j-6)*X3
        WVAL2= B(3*j-1) - D(9*j-5)*X1 - D(9*j-4)*X2 - D(9*j-3)*X3
        WVAL3= B(3*j  ) - D(9*j-2)*X1 - D(9*j-1)*X2 - D(9*j  )*X3
        do k= INL(j-1)+1, INL(j)
          i= IAL(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
          WVAL2= WVAL2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
          WVAL3= WVAL3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
        enddo
        do k= INU(j-1)+1, INU(j)
          i= IAU(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
          WVAL2= WVAL2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
          WVAL3= WVAL3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
        enddo

        WW(3*j-2,R)= WVAL1
        WW(3*j-1,R)= WVAL2
        WW(3*j  ,R)= WVAL3
      enddo
!$omp end parallel do 
      enddo
      enddo

      BNRM20= 0.d0
!$omp parallel do reduction(+:BNRM20)
      do i= 1, 3*N
        BNRM20= BNRM20+B(i)**2
      enddo
!$omp end parallel do 

      call MPI_Allreduce (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, SOLVER_COMM, ierr)
      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      ITER = 0
!C===

      call MPI_Barrier (SOLVER_COMM, ierr)
      S1_time= MPI_WTIME()
      do iter= 1, MAXIT
!      call MPI_BARRIER(SOLVER_COMM,ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===
        indexA= ZP
        indexB= Z
        indexC= R
!C
!C== Block SSOR

!$omp parallel do 
      do i= 1, 3*N
        WW(i,indexA)= WW(i,indexC)
        WW(i,indexB)= 0.d0
      enddo
!$omp end parallel do 

      do iterPRE= 1, iterPREmax

!C
!C-- FORWARD
        do iv= 1, NVECT
        do ip= 1,PEsmpTOT
          iS= STACKmc(ip-1,iv) + 1
          iE= STACKmc(ip  ,iv)
!$omp parallel do private(SW1,SW2,SW3,isL,ieL,j,k,X1,X2,X3)
        do i= iS, iE
          SW1= WW(3*i-2,indexA)
          SW2= WW(3*i-1,indexA)
          SW3= WW(3*i  ,indexA)
          isL= INL(i-1)+1
          ieL= INL(i)
          do j= isL, ieL
              k= IAL(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
            SW1= SW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            SW2= SW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            SW3= SW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
          enddo

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,indexA)= X1
          WW(3*i-1,indexA)= X2
          WW(3*i  ,indexA)= X3
        enddo
!$omp end parallel do 
        enddo
        enddo

        do iv= NVECT, 1, -1
        do ip= 1, PEsmpTOT
          iS= STACKmc(ip-1,iv) + 1
          iE= STACKmc(ip  ,iv)
!C
!C-- BACKWARD
!$omp parallel do private(isU,ieU,SW1,SW2,SW3,j,k,X1,X2,X3)
        do i= iS, iE
          isU= INU(i-1) + 1
          ieU= INU(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= isU, ieU
              k= IAU(j)
             X1= WW(3*k-2,indexA)
             X2= WW(3*k-1,indexA)
             X3= WW(3*k  ,indexA)
            SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
          enddo
          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,indexA)=  WW(3*i-2,indexA) - X1
          WW(3*i-1,indexA)=  WW(3*i-1,indexA) - X2
          WW(3*i  ,indexA)=  WW(3*i  ,indexA) - X3
        enddo
!$omp end parallel do 
        enddo
        enddo

!C
!C-- additive Schwartz
        if ( iterPRE .eq. 1 ) then
!$omp parallel do 
          do i= 1, N*3
            WW(i,indexB)= WW(i,indexA)
          enddo
!$omp end parallel do 
         else
!$omp parallel do 
          do i= 1, N*3
            WW(i,indexB)= WW(i,indexB) + WW(i,indexA)
          enddo
!$omp end parallel do 
        endif

        if (iterPRE.eq.iterPREmax) goto 750

        call SOLVER_SEND_RECV_3                                         &
     &     ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,        &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,indexB), 
     &       SOLVER_COMM,  my_rank)

      do jv= 1, NVECT
      do ip= 1,PEsmpTOT
        jS= STACKmc(ip-1,jv) + 1
        jE= STACKmc(ip  ,jv)
!$omp parallel do private(X1,X2,X3,WV1,WV2,WV3,k,i)
      do j= jS, jE
         X1= WW(3*j-2,indexB)
         X2= WW(3*j-1,indexB)
         X3= WW(3*j  ,indexB)
        WV1= WW(3*j-2,indexC) - D(9*j-8)*X1 - D(9*j-7)*X2 - D(9*j-6)*X3
        WV2= WW(3*j-1,indexC) - D(9*j-5)*X1 - D(9*j-4)*X2 - D(9*j-3)*X3
        WV3= WW(3*j  ,indexC) - D(9*j-2)*X1 - D(9*j-1)*X2 - D(9*j  )*X3
        do k= INL(j-1)+1, INL(j)
            i= IAL(k)
           X1= WW(3*i-2,indexB)
           X2= WW(3*i-1,indexB)
           X3= WW(3*i  ,indexB)
          WV1= WV1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
          WV2= WV2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
          WV3= WV3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
        enddo
        do k= INU(j-1)+1, INU(j)
            i= IAU(k)
           X1= WW(3*i-2,indexB)
           X2= WW(3*i-1,indexB)
           X3= WW(3*i  ,indexB)
          WV1= WV1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
          WV2= WV2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
          WV3= WV3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
        enddo
        WW(3*j-2,indexA)= WV1
        WW(3*j-1,indexA)= WV2
        WW(3*j  ,indexA)= WV3
      enddo
!$omp end parallel do 
      enddo
      enddo

 750  continue

      enddo
!C===
      
!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0

!$omp parallel do reduction(+:RHO0)
      do i= 1, 3*N
        RHO0= RHO0 + WW(i,R)*WW(i,Z)
      enddo
!$omp end parallel do 
      call MPI_Allreduce (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
     &                    MPI_SUM, SOLVER_COMM, ierr)
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then

!$omp parallel do 
        do i= 1, 3*N
          WW(i,P)= WW(i,Z)
        enddo
!$omp end parallel do 
       else
         BETA= RHO / RHO1

!$omp parallel do
        do i= 1, 3*N
          WW(i,P)= WW(i,Z) + BETA*WW(i,P)
        enddo
!$omp end parallel do 
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        

!C
!C-- INTERFACE data EXCHANGE
      call SOLVER_SEND_RECV_3                                           &
     &   ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,          &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,P) , SOLVER_COMM,     &
     &     my_rank)

!C
!C-- BEGIN calculation
      do jv= 1, NVECT
        jS= IVECT(jv-1) + 1
        jE= IVECT(jv)

      do ip= 1,PEsmpTOT
        jS= STACKmc(ip-1,jv) + 1
        jE= STACKmc(ip  ,jv)
!$omp parallel do private(X1,X2,X3,WVAL1,WVAL2,WVAL3,k,i)
      do j= jS, jE
           X1= WW(3*j-2,P)
           X2= WW(3*j-1,P)
           X3= WW(3*j  ,P)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,Q)= WVAL1
        WW(3*j-1,Q)= WVAL2
        WW(3*j  ,Q)= WVAL3

      enddo
!$omp end parallel do 
      enddo
      enddo

!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
!$omp parallel do reduction(+:C10)
      do i= 1, 3*N
        C10= C10 + WW(i,P)*WW(i,Q)
      enddo
!$omp end parallel do 

      call MPI_Allreduce (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
     &                    MPI_SUM, SOLVER_COMM, ierr)
      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===

!$omp parallel do 
      do i= 1, 3*N
         X(i  )= X (i  ) + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo
!$omp end parallel do 

      DNRM20= 0.d0
!$omp parallel do reduction(+:DNRM20)
      do i= 1, 3*N
        DNRM20= DNRM20 + WW(i,R)**2
      enddo
!$omp end parallel do 

      call MPI_Allreduce (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, SOLVER_COMM, ierr)
      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
#if 1
        if ((my_rank.eq.0) .and. (mod(iter,100)==0)) then
          write (*,1000) iter, RESID
        endif
#endif
 1000   format (i5, 1pe16.6)
! 1010   format (1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300

        RHO1 = RHO                                                             
      enddo
!C===

!C
!C-- INTERFACE data EXCHANGE
   30 continue
      E1_time= MPI_WTIME()

      COMPtime= (E1_time - S1_time)
      call MPI_Allreduce (COMPtime, CTsum, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_SUM, SOLVER_COMM, ierr)
      call MPI_Allreduce (COMPtime, CTmax, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_MAX, SOLVER_COMM, ierr)
      call MPI_Allreduce (COMPtime, CTmin, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_MIN, SOLVER_COMM, ierr)

      if (my_rank.eq.0) then
        write (*, 1000) ITER, RESID
        call MPI_COMM_SIZE (SOLVER_COMM, nPETOT, ierr )
        write (*,'(a,i8,3(1pe16.6))') 'elapsed', ITER, CTmax, CTmin,
     &                                           CTsum/dfloat(nPETOT)
        CTmax= CTmax/dfloat(ITER)
        CTmin= CTmin/dfloat(ITER)
        CTave= CTsum/dfloat(nPETOT)/dfloat(ITER)
        write (*,'(a,i8,3(1pe16.6))') 'iter   ', ITER, CTmax, CTmin, 
     &                                                 CTave
      endif

      call SOLVER_SEND_RECV_3                                           &
     &   ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,          &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

!$omp parallel do private(in)
      do i= 1, NP
        in= NtoO(i)
        WS(3*in-2)= X(3*i-2) 
        WS(3*in-1)= X(3*i-1) 
        WS(3*in  )= X(3*i  ) 
      enddo
!$omp end parallel do 

!$omp parallel do 
      do i= 1, NP
        X(3*i-2)= WS(3*i-2) 
        X(3*i-1)= WS(3*i-1) 
        X(3*i  )= WS(3*i  ) 
      enddo
!$omp end parallel do


      end subroutine        CG_3_SMP_novec
      end module     solver_CG_3_SMP_novec
