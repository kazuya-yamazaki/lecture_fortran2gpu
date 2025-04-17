!C
!C***
!C*** MAT_CON1
!C***
!C
      subroutine MAT_CON1
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer, dimension(:), allocatable :: IW3
      integer, dimension(:,:,:), allocatable :: IW2

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      allocate (IVECT   (0:2000,LEVELtot))
      allocate (OLDtoNEW(NP), NEWtoOLD(NP), IW(NP))
      allocate (NHYP(LEVELtot))

      IVECT   = 0
      OLDtoNEW= 0
      NEWtoOLD= 0
      IW      = 0
      NHYP    = 0

!C
!C +---------------+
!C | INIT_ORDERING |
!C +---------------+
!C===
      do lev= 1, LEVELtot
        NNlev= LEVELindex(lev) - LEVELindex(lev-1)

        if (NCOLORk.le.0) then
          call sRCM (NP, N, NU, NL, NHYP, INU, IAU, INL, IAL, IVECT,    &
     &               NEWtoOLD, OLDtoNEW, IW, lev, LEVELtot, LEVELindex, &
     &               my_rank)
         else
          if (lev.eq.1) then
            NCOLORx= NCOLORk
           else
            if (NNlev/NCOLORk.ge.10) then
              NCOLORx= NCOLORk
             else
              NCOLORx= NNlev/10 + 1
            endif
          endif
          call sMC  (NP, N, NU, NL, NHYP, INU, IAU, INL, IAL, IVECT,    &
     &               NEWtoOLD, OLDtoNEW, IW, lev, LEVELtot, LEVELindex, &
     &               my_rank, NCOLORx)
        endif
      enddo

!C
!C +-----------------+
!C | MATRIX transfer |
!C +-----------------+
!C===
      do i= N+1, NP
        OLDtoNEW(i)= i
        NEWtoOLD(i)= i
      enddo

      do j= 1, NL
        do i= 1, N
          IW(i) = IAL(NEWtoOLD(i),j)
        enddo
        do i= 1, N
          IAL(i,j) = IW(i)
        enddo
      enddo

      do j= 1, NU
        do i= 1, N
          IW(i) = IAU(NEWtoOLD(i),j)
        enddo
        do i= 1, N
          IAU(i,j) = IW(i)
        enddo
      enddo

      do i= 1, N
        IW(i) = INL(NEWtoOLD(i))
      enddo

      do i= 1, N
        INL(i) = IW(i)
      enddo

      do i= 1, N
        IW(i) = INU(NEWtoOLD(i))
      enddo

      do i= 1, N
        INU(i) = IW(i)
      enddo

      do j= 1, NL
        do i= 1, N
          if (IAL(i,j).eq.0) then
            IAL(i,j) = 0
           else
            IAL(i,j) = OLDtoNEW(IAL(i,j))
          endif
        enddo
      enddo

      do j= 1, NU
        do i= 1, N
          if ( IAU(i,j).eq.0) then
            IAU(i,j) = 0
           else
            IAU(i,j) = OLDtoNEW(IAU(i,j))
          endif
        enddo
      enddo

      do i= 1, NP
        IW= 0
        do k= 1, INL(i)
          IW(k)= IAL(i,k)
        enddo
        do k= 1, INU(i)
          IW(k+INL(i))= IAU(i,k)
        enddo
        icouK = INL(i) + INU(i)
        INL(i)= 0
        INU(i)= 0
        do k= 1, icouK
          jj= IW(k)
          if (jj.lt.i) then
            INL(i)= INL(i) + 1
            IAL(i,INL(i))= jj
           else
            INU(i)= INU(i) + 1
            IAU(i,INU(i))= jj
          endif
        enddo
      enddo
!C===

!C
!C +----------------------+
!C | CYCLIC multicoloring |
!C +----------------------+
!C===
      if (NCOLORk.ge.0) then
        NCOLORtot= NHYP(1)
       else
        NCOLORtot= -NCOLORk
      endif

      if (my_rank.eq.0) then
        write (*,'(a,i8)') 'color number: ', NCOLORtot 
      endif

      if (NCOLORk.ge.0) goto 1000
      
      allocate (ICHK(2000), OLDtoNEWmc(NP), NEWtoOLDmc(NP))
      allocate (IVECmc(0:2000), IVnew(NP))

      IW    = 0
      IVnew = 0
      NEWtoOLDmc= 0
      OLDtoNEWmc= 0

      do lev= 1, LEVELtot
      if (NHYP(lev).ne.0) then
        NCOLORtot= -NCOLORk
        if (lev.ge.2) NCOLORtot= NHYP(lev)/10 + 2

        IIs  = LEVELindex(lev-1)
        IImin= LEVELindex(lev-1)+1
        IImax= LEVELindex(lev)

  999   continue
        IVnew = 0
        IVECmc= 0
        do ic= 1, NHYP(lev)
          do ii= IVECT(ic-1,lev)+1, IVECT(ic,lev)
            IVnew(ii)= ic
          enddo
        enddo

        do ic= 1, NHYP(lev)
          ICHK(ic)= mod(ic-1,NCOLORtot) + 1
        enddo

        do i= IImin, IImax
          ic= IVnew(i)
          IVnew(i)= ICHK(ic)
        enddo
!C
!C-- CHECK dependency
        IFLAG= 0
        do i= IImin, IImax
          ic0= IVnew(i)
          do k= 1, INL(i)
            in= IAL(i,k)
            if (ic0.eq.IVnew(in)) IFLAG= 1
          enddo
          do k= 1, INU(i)
            in= IAU(i,k)
            if (ic0.eq.IVnew(in)) IFLAG= 1
          enddo

          if (IFLAG.eq.1) then
            NCOLORtot= NCOLORtot + 1
            if (my_rank.eq.0.and.lev.eq.1) then
              write (*,'(a,i8)') 'color number: ', NCOLORtot 
            endif
            goto 999
          endif
        enddo
!C
!C-- temp. ID
        icou= 0
        do ic= 1, NCOLORtot
          do i= IImin, IImax
            if (IVnew(i).eq.ic) then
              icou= icou + 1
              OLDtoNEWmc(i       )= icou + IIs
              NEWtoOLDmc(icou+IIs)= i
            endif
          enddo
          IVECT(ic,lev)= icou + IIs
        enddo
        NHYP(lev)= NCOLORtot
      endif      
      enddo

!C
!C-- matrix transfer
      do i= N+1, NP
        OLDtoNEWmc(i)= i
        NEWtoOLDmc(i)= i
      enddo

      do i= 1, N
        i1= OLDtoNEW  (i)
        i2= OLDtoNEWmc(i1)
        IVnew(i)= i2
      enddo

      do i= 1, N
        i1= IVnew(i)
        OLDtoNEW(i )= i1
        NEWtoOLD(i1)= i
      enddo

      do j= 1, NL
        do i= 1, N
          IW(i) = IAL(NEWtoOLDmc(i),j)
        enddo
        do i= 1, N
          IAL(i,j) = IW(i)
        enddo
      enddo

      do j= 1, NU
        do i= 1, N
          IW(i) = IAU(NEWtoOLDmc(i),j)
        enddo
        do i= 1, N
          IAU(i,j) = IW(i)
        enddo
      enddo

      do i= 1, N
        IW(i) = INL(NEWtoOLDmc(i))
      enddo

      do i= 1, N
        INL(i) = IW(i)
      enddo

      do i= 1, N
        IW(i) = INU(NEWtoOLDmc(i))
      enddo

      do i= 1, N
        INU(i) = IW(i)
      enddo

      do j= 1, NL
        do i= 1, N
          if (IAL(i,j).eq.0) then
            IAL(i,j) = 0
           else
            IAL(i,j) = OLDtoNEWmc(IAL(i,j))
          endif
        enddo
      enddo

      do j= 1, NU
        do i= 1, N
          if ( IAU(i,j).eq.0) then
            IAU(i,j) = 0
           else
            IAU(i,j) = OLDtoNEWmc(IAU(i,j))
          endif
        enddo
      enddo

      do i= 1, N
        IW= 0
        levi= LEVEL(NEWtoOLD(i))
        do k= 1, INL(i)
          IW(k)= IAL(i,k)
        enddo
        do k= 1, INU(i)
          IW(k+INL(i))= IAU(i,k)
        enddo
        icouK = INL(i) + INU(i)
        INL(i)= 0
        INU(i)= 0
        do k= 1, icouK
          jj= IW(k)
          levj= LEVEL(NEWtoOLD(jj))          

          if (jj.gt.N) then
            if (levj.ge.levi) then
              INU(i)= INU(i) + 1
              IAU(i,INU(i))= jj
             else
              INL(i)= INL(i) + 1
              IAL(i,INL(i))= jj
            endif
           else    
            if (jj.lt.i) then
              INL(i)= INL(i) + 1
              IAL(i,INL(i))= jj
             else
              INU(i)= INU(i) + 1
              IAU(i,INU(i))= jj
            endif
          endif

        enddo
      enddo

 1000 continue
 2000 continue
!C===

!C
!C +---------+
!C | for SMP |
!C +---------+
!C===
      if (my_rank.eq.0) write (*,'(/a)') 
      allocate (STACKmcG(0:PEsmpTOT))

      NHYPmax= -100
      do lev= 1, LEVELtot
        NHYPmax= max(NHYPmax, NHYP(lev))
      enddo
      allocate (STACKmc(0:PEsmpTOT, NHYPmax, LEVELtot))
      STACKmc= 0
      do lev= 1, LEVELtot
        do ic= 1, NHYP(lev)
          NN0= IVECT(ic,lev) - IVECT(ic-1,lev)
          NN= NN0/PEsmpTOT
          NR= NN0 - NN*PEsmpTOT
          STACKmc(0,ic,lev)= IVECT(ic-1,lev)
          do ip= 1, PEsmpTOT
            STACKmc(ip,ic,lev)= NN
            if (ip.le.NR) STACKmc(ip,ic,lev)= NN + 1
          enddo
          do ip= 1, PEsmpTOT
            STACKmc(ip,ic,lev)= STACKmc(ip  ,ic,lev) + 
     &                          STACKmc(ip-1,ic,lev)
          enddo
        enddo
      enddo
!C===


!C
!C +-------------------+
!C | FINAL RE-ORDERING |
!C +-------------------+
!C===
      NLmax= -NP
      NUmax= -NP

      do i= 1, N
        NLmax= max(NLmax, INL(i))
        NUmax= max(NUmax, INU(i))
      enddo

      deallocate (OLDtoNEWmc, NEWtoOLDmc)
      allocate (OLDtoNEWmc(0:NP), NEWtoOLDmc(0:NP))
      NEWtoOLDmc= 0
      OLDtoNEWmc= 0

      do i= N+1, NP
        OLDtoNEWmc(i)= i
        NEWtoOLDmc(i)= i
      enddo

      allocate (IW3(0:PEsmpTOT))
      STACKmcG= 0
      IW3= 0

      IW= -1
      do lev= 1, LEVELtot
        do ic= 1, NHYP(lev)
          do ip= 1, PEsmpTOT
          do i = STACKmc(ip-1,ic,lev)+1, STACKmc(ip,ic,lev)
            IW(i)= ip
            STACKmcG(ip)= STACKmcG(ip) + 1
          enddo
          enddo
        enddo
      enddo

      do ip= 1, PEsmpTOT
        STACKmcG(ip)= STACKmcG(ip-1) + STACKmcG(ip)
      enddo

      do lev= 1, LEVELtot
        do ic= 1, NHYP(lev)
          do ip= 1, PEsmpTOT
          do i = STACKmc(ip-1,ic,lev)+1, STACKmc(ip,ic,lev)
            IW3(ip)= IW3(ip) + 1
            icou= STACKmcG(ip-1) + IW3(ip)
            OLDtoNEWmc(i   )= icou
            NEWtoOLDmc(icou)= i
          enddo
          enddo
        enddo
      enddo

      allocate (IW2(PEsmpTOT, 0:NHYPmax, LEVELtot))
      IW2= 0
      do ip= 1, PEsmpTOT
        IW2(ip,0,1)= STACKmcG(ip-1)
      enddo

      do ip= 1, PEsmpTOT
        do lev= 1, LEVELtot
        do ic= 1, NHYP(lev)
            istart= STACKmc(ip-1,ic,lev)
            iend  = STACKmc(ip  ,ic,lev)
            IW2(ip,ic,lev)= IW2(ip,ic-1,lev) + iend - istart 
        enddo
        if (lev.lt.LEVELtot) then
          IW2(ip,0,lev+1)= IW2(ip,NHYP(lev),lev)
        endif
        enddo
      enddo

      deallocate (STACKmc)
      allocate (STACKmc(PEsmpTOT, 0:NHYPmax, LEVELtot))
      STACKmc= IW2

      do i= 1, N
        in1= OLDtoNEW  (i)
        in2= OLDtoNEWmc(in1)
        OLDtoNEW(i  )= in2
        NEWtoOLD(in2)= i
      enddo

      do i= 1, N
        IW(i)= INL(i)
      enddo
      do i= 1, N
        INL(i)= IW(NEWtoOLDmc(i))
      enddo
      do i= 1, N
        IW(i)= INU(i)
      enddo
      do i= 1, N
        INU(i)= IW(NEWtoOLDmc(i))
      enddo

      do j= 1, NLmax
        do i= 1, N
          IW(i)= IAL(i,j)
        enddo
        do i= 1, N
          IAL(i,j)= OLDtoNEWmc(IW(NEWtoOLDmc(i)))
        enddo
      enddo

      do j= 1, NUmax
        do i= 1, N
          IW(i)= IAU(i,j)
        enddo
        do i= 1, N
          IAU(i,j)= OLDtoNEWmc(IW(NEWtoOLDmc(i)))
        enddo
      enddo
      deallocate (OLDtoNEWmc, NEWtoOLDmc)
!C===

!C
!C +-------------------+
!C | compressed MATRIX |
!C +-------------------+
!C===
      allocate (indexL(0:NP), indexU(0:NP))
      indexL(0)= 0
      indexU(0)= 0

      if (FTflag.eq.1) then
        do lev= 1, LEVELtot
          do ic= 1, NHYP(lev)
!$omp parallel do private(ip,i)
            do ip= 1, PEsmpTOT
            do i = STACKmc(ip,ic-1,lev)+1, STACKmc(ip,ic,lev)
              indexL(i)= 0
              indexU(i)= 0
            enddo
            enddo
          enddo
        enddo
       else
        indexL= 0
        indexU= 0
      endif

      do i= 1, N
        indexL(i)= indexL(i-1) + INL(i)
        indexU(i)= indexU(i-1) + INU(i)
      enddo

      NPL= indexL(N)
      NPU= indexU(N)

      allocate (itemL(NPL), itemU(NPU))

      if (FTflag.eq.1) then
        do lev= 1, LEVELtot
          do ic= 1, NHYP(lev)
!$omp parallel do private(ip,i,j,isL,ieL,isU,ieU)
            do ip= 1, PEsmpTOT
            do i = STACKmc(ip,ic-1,lev)+1, STACKmc(ip,ic,lev)
              isL= indexL(i-1)+1
              ieL= indexL(i)
              do j= isL, ieL
                itemL(j)= 0
              enddo

              isU= indexU(i-1)+1
              ieU= indexU(i)
              do j= isU, ieU
                itemU(j)= 0
              enddo
            enddo
            enddo
          enddo
        enddo
       else
        itemL= 0
        itemU= 0
      endif

      do i= 1, N
        do k= 1, INL(i)
                kk = k + indexL(i-1)
          itemL(kk)=        IAL(i,k)
        enddo
        do k= 1, INU(i)
                kk = k + indexU(i-1)
          itemU(kk)=        IAU(i,k)
        enddo
      enddo

      deallocate (INL, INU, IAL, IAU)

      do neib= 1, NEIBPEtot
        do k= EXPORT_index(neib-1)+1, EXPORT_index(neib)
          kk= EXPORT_ITEM(k)
          EXPORT_ITEM(k)= OLDtoNEW(kk)
        enddo
      enddo
!C===
      end subroutine MAT_CON1
