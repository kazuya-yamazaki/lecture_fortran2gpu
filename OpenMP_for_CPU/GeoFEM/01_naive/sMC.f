!C    
!C***
!C*** sMC
!C***
!C
!C    Multicolor Ordering Method 
!C    

      subroutine sMC (NP, N, NL, NU, NCOLORtot, INL, IAL, INU, IAU,     &
     &                COLORindex, NEWtoOLD, OLDtoNEW)

      implicit REAL*8(A-H,O-Z)

      integer, dimension(NP)   :: INL, INU, NEWtoOLD, OLDtoNEW
      integer, dimension(0:NP) :: COLORindex
      integer, dimension(NP,NL):: IAL
      integer, dimension(NP,NU):: IAU

      integer, dimension(:)  , allocatable :: IW, INLw, INUw
      integer, dimension(:,:), allocatable :: IALw, IAUw

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      allocate (IW(NP))
      IW= 0

      NCOLORk  = NCOLORtot

       INmin= N
      NODmin= 0
      do i= 1, N
        icon= INL(i) + INU(i)
        if (icon.lt.INmin) then
           INmin= icon
          NODmin= i
        endif
      enddo

      OLDtoNEW(NODmin)= 1
      NEWtoOLD(     1)= NODmin

      IW        = 0
      IW(NODmin)= 1

      ITEMcou= N/NCOLORk
!C===

!C
!C +---------------+
!C | MULTIcoloring |
!C +---------------+
!C===
      icou = 1
      icouK= 1

      do icol= 1, N
        NCOLORk= icol
        do i= 1, N
          if (IW(i).le.0) IW(i)= 0
        enddo
        do i= 1, N
!C
!C-- already COLORED
          if (IW(i).eq.icol) then
            do k= 1, INL(i)
              ik= IAL(i,k)
              if (IW(ik).le.0) IW(ik)= -1
            enddo
            do k= 1, INU(i)
              ik= IAU(i,k)
              if (IW(ik).le.0) IW(ik)= -1
            enddo
          endif
!C
!C-- not COLORED
          if (IW(i).eq.0) then
            icou = icou  + 1
            icouK= icouK + 1
            IW(i)= icol
            do k= 1, INL(i)
              ik= IAL(i,k)
              if (IW(ik).le.0) IW(ik)= -1
            enddo
            do k= 1, INU(i)
              ik= IAU(i,k)
              if (IW(ik).le.0) IW(ik)= -1
            enddo
          endif
          if (icou .eq.N)       goto 200
          if (icouK.eq.ITEMcou) goto 100
        enddo      
 100    continue
        icouK= 0
      enddo

 200  continue
!C===

!C
!C +----------------+
!C | FINAL COLORING |
!C +----------------+
!C===
      NCOLORtot= NCOLORk
      COLORindex= 0
      icoug= 0
      do ic= 1, NCOLORtot
        icou= 0
        do i= 1, N
          if (IW(i).eq.ic) then
            icou = icou + 1
            icoug= icoug + 1
            NEWtoOLD(icoug)= i
            OLDtoNEW(i    )= icoug
          endif
        enddo
        COLORindex(ic)= icou
      enddo

      do ic= 1, NCOLORtot
        COLORindex(ic)= COLORindex(ic-1) + COLORindex(ic)
      enddo
!C===

!C
!C +-----------------+
!C | MATRIX transfer |
!C +-----------------+
!C===
      allocate (INLw(NP), INUw(NP), IALw(NP,NL), IAUw(NP,NU))
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

      do i= 1, NP
        IW(i) = INL(NEWtoOLD(i))
      enddo

      do i= 1, NP
        INLw(i) = IW(i)
      enddo

      do i= 1, NP
        IW(i) = INU(NEWtoOLD(i))
      enddo

      do i= 1, NP
        INUw(i) = IW(i)
      enddo

      do j= 1, NL
        do i= 1, NP
          if (IAL(i,j).eq.0) then
            IALw(i,j) = 0
           else
            IALw(i,j) = OLDtoNEW(IAL(i,j))
          endif
        enddo
      enddo

      do j= 1, NU
        do i= 1, NP
          if ( IAU(i,j).eq.0) then
            IAUw(i,j) = 0
           else
            IAUw(i,j) = OLDtoNEW(IAU(i,j))
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
!C===
      deallocate (IW, INLw, INUw, IALw, IAUw)

      return
      end



