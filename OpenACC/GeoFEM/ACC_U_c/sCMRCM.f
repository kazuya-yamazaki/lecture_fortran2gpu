!C    
!C***
!C*** sCM-RCM
!C***
!C
!C    Cyclic Multicolored RCM
!C    

      subroutine sCMRCM (NP, N, NL, NU, NHYP, INL, IAL, INU, IAU, 
     &                   COLORindex, NEWtoOLD, OLDtoNEW, my_rank)

      implicit REAL*8(A-H,O-Z)

      integer :: CMtot
      integer, dimension(NP)   :: INL, INU, NEWtoOLD, OLDtoNEW
      integer, dimension(0:NP) :: COLORindex
      integer, dimension(NP,NL):: IAL
      integer, dimension(NP,NU):: IAU

      integer, dimension(:), allocatable :: IW

      integer, dimension(:)  , allocatable :: INLw, INUw
      integer, dimension(:,:), allocatable :: IALw, IAUw

      integer, dimension(:)  , allocatable :: ICHK, IVECcm, IVnew
      integer, dimension(:)  , allocatable :: OLDtoNEWcm, NEWtoOLDcm

!C
!C-- INIT.
      NCOLORtot= -NHYP

      allocate (IW(NP))
      allocate (INLw(NP), INUw(NP), IALw(NP,NL), IAUw(NP,NU))

      IW   = 0
      COLORindex= 0

      COLORindex(0)= 0
      COLORindex(1)= 0

      NHYP= 0
      KCT = 0

!C
!C-- RCM ordering
      do i= 1, N
        if (INL(i).eq.0) then
          COLORindex(1) = COLORindex(1) + 1
          KCT           = KCT + 1
          NEWtoOLD(KCT) = i
          IW (i)        = -1
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

      NHYP = NHYP + 1
      KCT0 = KCT
      do i = KMIN, KMAX
        if (IW(i).eq.INL(i)) then
              KCT = KCT + 1
          IW (i)  = -1
          NEWtoOLD(KCT)= I
        endif
      enddo

      KC                 = KCT - KCT0
      COLORindex(NHYP+1) = KCT

      if (KC.ne.0 .and. NHYP.le.N) goto 20

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

!C
!C +-------------------+
!C | Cyclic-ReORDERING |
!C +-------------------+
!C===
      if (NCOLORtot.gt.NHYP) then
        NCOLORtot= NHYP
        goto 2001
      endif

 999  continue
      if (my_rank.eq.0) write (*,'(a,i8)') 'color number: ', NCOLORtot

      if (allocated(ICHK      )) deallocate (ICHK)
      if (allocated(OLDtoNEWcm)) deallocate (OLDtoNEWcm)
      if (allocated(NEWtoOLDcm)) deallocate (NEWtoOLDcm)
      if (allocated(IVECcm    )) deallocate (IVECcm    )
      if (allocated(IVnew     )) deallocate (IVnew     )

      allocate (ICHK(NHYP), OLDtoNEWcm(NP), NEWtoOLDcm(NP))
      allocate (IVECcm(0:NCOLORtot), IVnew(NP))

      IW    = 0
      IVECcm= 0
      IVnew = 0

      do i= 1, NP
        OLDtoNEWcm(i)= i
        NEWtoOLDcm(i)= i
      enddo

      do i= 1, NCOLORtot
        ICHK = 0
        icou = i - NCOLORtot
        icolt= 0
        do k= 1, NHYP
          icou= icou + NCOLORtot
          if (icou.gt.NHYP) exit
          IW(i)= IW(i) + COLORindex(icou)-COLORindex(icou-1)
          icolt= icolt + 1
          ICHK(icolt)= icou
        enddo

        CMtot= icolt
        ICHK1= 0
        icoug= 0

        do ic1= 1, CMtot
          ic2= ICHK(ic1)
          do k= COLORindex(ic2-1)+1, COLORindex(ic2)
            icoug    = icoug + 1
            IVnew(icoug+IVECcm(i-1))= k
          enddo
        enddo
        IVECcm(i)= IVECcm(i-1) + icoug
      enddo

!C
!C-- CHECK dependency
      IFLAG= 0
      do i= 1, NCOLORtot
        iS= IVECcm(i-1) + 1
        iE= IVECcm(i)
        IW= 0
        do j= iS, iE
          in= IVnew(j)
          IW(in)= 1
        enddo

        do j= iS, iE
          in= IVnew(j)
          do k= 1, INL(in)
            ip= IAL(in,k)
            if (ip.ne.0) then
              if (IW(ip).eq.1) IFLAG= 1
            endif
          enddo
          do k= 1, INU(in)
            ip= IAU(in,k)
            if (ip.ne.0) then
              if (IW(ip).eq.1) IFLAG= 1
            endif           
          enddo
        enddo
      enddo

      if (IFLAG.eq.1) then
        NCOLORtot= NCOLORtot + 1
        deallocate (ICHK, OLDtoNEWcm, NEWtoOLDcm, IVECcm, IVnew)
        goto 999
      endif

 2000  continue
      do i= 1, N
        in= IVnew(i)
        OLDtoNEWcm(in)= i
      enddo

      do i= 1, N
        in= OLDtoNEWcm(i)
        NEWtoOLDcm(in)= i
      enddo
!C===

!C
!C +-----------------+
!C | MATRIX transfer |
!C +-----------------+
!C===
      INLw= 0
      INUw= 0
      IALw= 0
      IAUw= 0

      do j= 1, NL
        do i= 1, NP
          IW(i) = IAL(NEWtoOLDcm(i),j)
        enddo
        do i= 1, NP
          IAL(i,j) = IW(i)
        enddo
      enddo

      do j= 1, NU
        do i= 1, NP
          IW(i) = IAU(NEWtoOLDcm(i),j)
        enddo
        do i= 1, NP
          IAU(i,j) = IW(i)
        enddo
      enddo

      do i= 1, NP
        IW(i) = INL(NEWtoOLDcm(i))
      enddo

      do i= 1, NP
        INLw(i) = IW(i)
      enddo

      do i= 1, NP
        IW(i) = INU(NEWtoOLDcm(i))
      enddo

      do i= 1, NP
        INUw(i) = IW(i)
      enddo

      do j= 1, NL
        do i= 1, NP
          if (IAL(i,j).eq.0) then
            IALw(i,j) = 0
           else
            IALw(i,j) = OLDtoNEWcm(IAL(i,j))
          endif
        enddo
      enddo

      do j= 1, NU
        do i= 1, NP
          if ( IAU(i,j).eq.0) then
            IAUw(i,j) = 0
           else
            IAUw(i,j) = OLDtoNEWcm(IAU(i,j))
          endif
        enddo
      enddo

      INL= 0
      INU= 0
      IAL= 0
      IAU= 0

      do i= 1, NP
        jL= 0
        jU= 0
        do j= 1, INLw(i)
          if (IALw(i,j).gt.i) then
            jU= jU + 1
            IAU(i,jU)= IALw(i,j)
           else
            jL= jL + 1
            IAL(i,jL)= IALw(i,j)
          endif
         enddo

        do j= 1, INUw(i)
          if (IAUw(i,j).gt.i) then
            jU= jU + 1
            IAU(i,jU)= IAUw(i,j)
           else
            jL= jL + 1
            IAL(i,jL)= IAUw(i,j)
          endif
         enddo

         INL(i)= jL
         INU(i)= jU
      enddo

      do i= 1, N
        ic1= OLDtoNEW(i)
        ic2= OLDtoNEWcm(ic1)
        OLDtoNEW(i)= ic2
      enddo

      do i= 1, N
        in= OLDtoNEW(i)
        NEWtoOLD(in)= i
      enddo

      COLORindex= 0
      do icol= 1, NCOLORtot
        COLORindex(icol)= IVECcm(icol)
      enddo
!C===
      deallocate (ICHK, OLDtoNEWcm, NEWtoOLDcm, IVECcm, IVnew)
      NHYP= NCOLORtot

 2001 continue

      deallocate (IW, INLw, INUw, IALw, IAUw)

      return
      end







