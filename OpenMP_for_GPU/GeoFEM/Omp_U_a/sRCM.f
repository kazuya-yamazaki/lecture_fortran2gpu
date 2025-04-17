!C    
!C***
!C*** sRCM
!C***
!C
!C    Reverse Cuthill-McKee
!C    

      subroutine sRCM (NP, N, NL, NU, NO, INL, IAL, INU, IAU, IVECT,    &
     &                 NEWtoOLD, OLDtoNEW, my_rank)

      implicit REAL*8(A-H,O-Z)

      integer, dimension(NP)   :: INL, INU, NEWtoOLD, OLDtoNEW
      integer, dimension(0:NP) :: IVECT
      integer, dimension(NP,NL):: IAL
      integer, dimension(NP,NU):: IAU

      integer, dimension(:), allocatable :: IW

!C
!C-- INIT.
      allocate (IW(NP))
      IW   = 0
      IVECT= 0

      IVECT(0)= 0
      IVECT(1)= 0
       NO   = 0
      KCT   = 0

!C
!C-- RCM ordering
      do i= 1, N
        if (INL(i).eq.0) then
          IVECT(1) = IVECT(1) + 1
          KCT      = KCT + 1
          NEWtoOLD(KCT) = i
          IW (i)   = -1
         else
          IW(i) = 0
       endif
      enddo  

      KC  = KCT
      KCT0= 0

   20 continue
      KMIN = N
      KMAX = 1
      do j= 1, KC
        JC = NEWtoOLD(KCT0+j)
        JN = INU(JC)
        do i= 1, JN
             II = IAU(JC,i)
          if (II.le.N) then
            IW(II)= IW(II) + 1
             KMIN = min0(II,KMIN)
             KMAX = max0(II,KMAX)
          endif
        enddo
      enddo

        NO = NO + 1
      KCT0 = KCT
      do i = KMIN, KMAX
        if (IW(i).eq.INL(i)) then
              KCT = KCT + 1
          IW (i)  = -1
          NEWtoOLD(KCT)= I
        endif
      enddo

       KC         = KCT - KCT0
      IVECT(NO+1) = KCT

      if (KC.ne.0 .and. NO.le.N) goto 20

!C
!C-- NEWtoOLD, OLDtoNEW array
      do i= 1, N
        i0= NEWtoOLD(i)
        OLDtoNEW(i0)= i
      enddo

      do j= 1, NL
        do i= 1, NP
          IW(i) = IAL(NEWtoOLD(i),j)
        enddo
        do i= 1, NP
          IAL(i,j) = IW(i)
        enddo
      enddo

      do j= 1, NU
        do i= 1, NP
          IW(i) = IAU(NEWtoOLD(i),j)
        enddo
        do i= 1, NP
          IAU(i,j) = IW(i)
        enddo
      enddo

!C
!C-- MATRIX transfer
      do i= 1, NP
        IW(i) = INL(NEWtoOLD(i))
      enddo

      do i= 1, NP
        INL(i) = IW(i)
      enddo

      do i= 1, NP
        IW(i) = INU(NEWtoOLD(i))
      enddo

      do i= 1, NP
        INU(i) = IW(i)
      enddo

      do j= 1, NL
        do i= 1, NP
          if (IAL(i,j).eq.0) then
            IAL(i,j) = 0
           else
            IAL(i,j) = OLDtoNEW(IAL(i,j))
          endif
        enddo
      enddo

      do j= 1, NU
        do i= 1, NP
          if ( IAU(i,j).eq.0) then
            IAU(i,j) = 0
           else
            IAU(i,j) = OLDtoNEW(IAU(i,j))
          endif
        enddo
      enddo

      return
      end







