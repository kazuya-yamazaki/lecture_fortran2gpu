!C
!C***
!C*** MAT_CON1
!C***
!C
      subroutine MAT_CON1
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer, dimension(0:1024) :: IW3
      integer, dimension(:,:), allocatable :: IW2

      call INIT_ORDERING
      call PREP_SMP
      call FINAL_ORDERING
      call ONEdim_ORDERING

      if (PETOT.ne.1) then
        do k= 1, NOD_STACK_EXPORT(NEIBPETOT)
          i= NOD_EXPORT(k)
          in= OLDtoNEW(i)
          NOD_EXPORT_NEW(k)= in
        enddo

        do k= 1, NOD_STACK_EXPORT(NEIBPETOT)
          NOD_EXPORT(k)= NOD_EXPORT_NEW(k)
        enddo
      endif

      contains
!C
!C***
!C*** INIT_ORDERING
!C***
!C
      subroutine INIT_ORDERING

      if (allocated(IW)) deallocate (IW)
      allocate (IW(NP),IVECT(0:NP))
      allocate (OLDtoNEW(NP), NEWtoOLD(NP))

      IVECT   = 0
      IW      = 0

      do i= 1, NP
        OLDtoNEW(i)= i
        NEWtoOLD(i)= i
      enddo

      if (NCOLORk.eq.0.or.NCOLORk.eq.-1) then
        if (my_rank.eq.0) write (*,'(a)') '### RCM'
        call sRCM  (NP, N, NL, NU, NHYP, INL, IAL, INU, IAU, IVECT,     &
     &              NEWtoOLD, OLDtoNEW, my_rank)
      endif
      if (NCOLORk.lt.-1) then
        if (my_rank.eq.0) write (*,'(a)') '### CM-RCM'
        NHYP= NCOLORk
        call sCMRCM  (NP, N, NL, NU, NHYP, INL, IAL, INU, IAU, IVECT,   &
     &                NEWtoOLD, OLDtoNEW, my_rank)  
      endif
      if (NCOLORk.gt.1) then
        if (my_rank.eq.0) write (*,'(a)') '### Multicoloring'
        NHYP= NCOLORk
        call sMC  (NP, N, NL, NU, NHYP, INL, IAL, INU, IAU, IVECT,      &
     &             NEWtoOLD, OLDtoNEW)
      endif

      NLmax= -NP
      NUmax= -NP

      do i= 1, NP
        NLmax= max(NLmax, INL(i))
        NUmax= max(NUmax, INU(i))
      enddo

      NCOLORtot= NHYP
      if (my_rank.eq.0) write (*,'(a,i8)') 'color number: ', NCOLORtot

      end subroutine INIT_ORDERING
!C
!C***
!C*** PREP_SMP
!C***
!C
      subroutine PREP_SMP
!C
!C-- SMP stack
      allocate (STACKmcG(0:PEsmpTOT))
      STACKmcG= 0

      icon= N/PEsmpTOT
      ir  = N - icon*PEsmpTOT
      do ip= 1, PEsmpTOT
        STACKmcG(ip)= icon
      enddo
      do ip= 1, ir
        STACKmcG(ip)= icon + 1
      enddo

      do ip= 1, PEsmpTOT
        STACKmcG(ip)= STACKmcG(ip-1) + STACKmcG(ip)
      enddo

!C
!C-- SMP stack
      if (allocated(STACKmc)) deallocate (STACKmc)
      allocate (STACKmc(0:PEsmpTOT, NCOLORtot))
      STACKmc = 0

      do ic= 1, NCOLORtot
        iS= IVECT(ic-1)
        iE= IVECT(ic  )
        icon= (iE-iS)/PEsmpTOT
        ir  = (iE-iS)-icon*PEsmpTOT

        do ip= 1, PEsmpTOT
          STACKmc(ip,ic)= icon
        enddo
        do ip= 1, ir
          STACKmc(ip,ic)= icon + 1
        enddo
      enddo

      do ic= 1, NCOLORtot
        if (ic.ne.1) STACKmc(0,ic)= STACKmc(PEsmpTOT,ic-1)
      do ip= 1, PEsmpTOT
        STACKmc(ip,ic)= STACKmc(ip-1,ic) + STACKmc(ip,ic)
      enddo
      enddo

      end subroutine PREP_SMP

!C
!C***
!C*** FINAL_ORDERING
!C***
!C
      subroutine FINAL_ORDERING
      end subroutine FINAL_ORDERING

!C
!C***
!C*** ONEdim_ORDERING
!C***
!C
      subroutine ONEdim_ORDERING
      
      if (allocated(indexL)) deallocate (indexL)
      if (allocated(indexU)) deallocate (indexU)
      allocate (indexL(0:NP), indexU(0:NP))      

      indexL(0)= 0
      indexU(0)= 0

      if (FTflag.eq.1) then
          do jv= 1, NCOLORtot
            jS= IVECT(jv-1) + 1
            jE= IVECT(jv)
!$omp parallel do private(ip,j, jS, jE)
            do ip= 1, PEsmpTOT
              jS= STACKmc(ip-1,jv) + 1
              jE= STACKmc(ip,  jv)
            do j= jS, jE
              indexL(j)= 0
              indexU(j)= 0
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
      if (allocated(itemL)) deallocate (itemL)
      if (allocated(itemU)) deallocate (itemU)
      allocate (itemL(NPL), itemU(NPU))

      if (FTflag.eq.1) then
          do jv= 1, NCOLORtot
            jS= IVECT(jv-1) + 1
            jE= IVECT(jv)
!$omp parallel do private(jS,jE,j,jsL,jeL,jsU,jeU,k)
            do ip= 1, PEsmpTOT
              jS= STACKmc(ip-1,jv) + 1
              jE= STACKmc(ip  ,jv)
            do j = jS, jE
              jsL= indexL(j-1)+1
              jeL= indexL(j)
              do k= jsL, jeL
                itemL(k)= 0
              enddo

              jsU= indexU(j-1)+1
              jeU= indexU(j)
              do k= jsU, jeU
                itemU(k)= 0
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
      end subroutine ONEdim_ORDERING

        subroutine FIND_TS_NODE (ip1,ip2)
        if (ip1.gt.ip2) then
          do kk= 1, INL(ip1)
            if (ip2.eq.IAL(ip1,kk)) return
          enddo
          icou= INL(ip1) + 1
          IAL(ip1,icou)= ip2
          INL(ip1     )= icou
          return
        endif

        if (ip2.gt.ip1) then
          do kk= 1, INU(ip1)
            if (ip2.eq.IAU(ip1,kk)) return
          enddo
          icou= INU(ip1) + 1
          IAU(ip1,icou)= ip2
          INU(ip1     )= icou
          return
        endif

        end subroutine FIND_TS_NODE

      end subroutine MAT_CON1
