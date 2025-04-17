!C
!C***
!C*** MAT_ASS_MAIN
!C***
!C
      subroutine MAT_ASS_MAIN 
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(  8) :: nodLOCAL


      call MPI_Barrier (SOLVER_COMM, ierr)
      S1_time= MPI_WTIME()
!C
!C +------------------+
!C | ELEMENT Coloring |
!C +------------------+
!C===
      allocate (ELMCOLORindex(0:NP))
      allocate (ELMCOLORitem (ICELTOT))
      if (allocated (IWKX)) deallocate (IWKX)
      allocate (IWKX(0:NP,3))

      IWKX= 0
      icou= 0
      do icol= 1, NP
        do i= 1, NP
          IWKX(i,1)= 0
        enddo
        do icel= 1, ICELTOT
          if (IWKX(icel,2).eq.0) then
            in1= ICELNOD(icel,1)
            in2= ICELNOD(icel,2)
            in3= ICELNOD(icel,3)
            in4= ICELNOD(icel,4)
            in5= ICELNOD(icel,5)
            in6= ICELNOD(icel,6)
            in7= ICELNOD(icel,7)
            in8= ICELNOD(icel,8)

            ip1= IWKX(in1,1)
            ip2= IWKX(in2,1)
            ip3= IWKX(in3,1)
            ip4= IWKX(in4,1)
            ip5= IWKX(in5,1)
            ip6= IWKX(in6,1)
            ip7= IWKX(in7,1)
            ip8= IWKX(in8,1)

            isum= ip1 + ip2 + ip3 + ip4 + ip5 + ip6 + ip7 + ip8
            if (isum.eq.0) then 
              icou= icou + 1
              IWKX(icol,3)= icou
              IWKX(icel,2)= icol
              ELMCOLORitem(icou)= icel

              IWKX(in1,1)= 1
              IWKX(in2,1)= 1
              IWKX(in3,1)= 1
              IWKX(in4,1)= 1
              IWKX(in5,1)= 1
              IWKX(in6,1)= 1
              IWKX(in7,1)= 1
              IWKX(in8,1)= 1
              if (icou.eq.ICELTOT) goto 100            
            endif
          endif
        enddo
      enddo

 100  continue
      ELMCOLORtot= icol
      IWKX(0          ,3)= 0
      IWKX(ELMCOLORtot,3)= ICELTOT

      do icol= 0, ELMCOLORtot
        ELMCOLORindex(icol)= IWKX(icol,3)
      enddo

!      write (*,'(a,2i8)') '### Number of Element Colors', 
!     &                     my_rank, ELMCOLORtot
      deallocate (IWKX)
      X1_time= MPI_WTIME()
!C===

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00

      if (allocated(AL)) deallocate (AL)
      if (allocated(AU)) deallocate (AU)
      if (allocated(D )) deallocate (D )
      if (allocated(B )) deallocate (B )
      if (allocated(X )) deallocate (X )
      if (allocated(WW)) deallocate (WW)
      if (allocated(ALUG)) deallocate (ALUG)
      if (allocated(WS)) deallocate (WR)

      allocate (AL(9*NPL), AU(9*NPU), ALUG(9*NP))
      allocate (D(9*NP), B(3*NP), X(3*NP), WW(3*NP,4))
      allocate (WS(3*NP), WR(3*NP))

      if (FTflag.eq.1) then
          do jv= 1, NHYP
            jS= IVECT(jv-1) + 1
            jE= IVECT(jv)
!$omp parallel do private(ip,j, jS, jE)
            do ip= 1, PEsmpTOT
              jS= STACKmc(ip-1,jv) + 1
              jE= STACKmc(ip  ,jv)
            do j= jS, jE
              X(3*j-2)= 0.d0
              X(3*j-1)= 0.d0
              X(3*j  )= 0.d0
              B(3*j-2)= 0.d0
              B(3*j-1)= 0.d0
              B(3*j  )= 0.d0
              D(9*j-8)= 0.d0
              D(9*j-7)= 0.d0
              D(9*j-6)= 0.d0
              D(9*j-5)= 0.d0
              D(9*j-4)= 0.d0
              D(9*j-3)= 0.d0
              D(9*j-2)= 0.d0
              D(9*j-1)= 0.d0
              D(9*j  )= 0.d0
              ALUG(9*j  )= 0.d0
              ALUG(9*j-8)= 0.d0
              ALUG(9*j-7)= 0.d0
              ALUG(9*j-6)= 0.d0
              ALUG(9*j-5)= 0.d0
              ALUG(9*j-4)= 0.d0
              ALUG(9*j-3)= 0.d0
              ALUG(9*j-2)= 0.d0
              ALUG(9*j-1)= 0.d0
              ALUG(9*j  )= 0.d0
              WW(3*j-2,1)= 0.d0
              WW(3*j-1,1)= 0.d0
              WW(3*j  ,1)= 0.d0
              WW(3*j-2,2)= 0.d0
              WW(3*j-1,2)= 0.d0
              WW(3*j  ,2)= 0.d0
              WW(3*j-2,3)= 0.d0
              WW(3*j-1,3)= 0.d0
              WW(3*j  ,3)= 0.d0
              WW(3*j-2,4)= 0.d0
              WW(3*j-1,4)= 0.d0
              WW(3*j  ,4)= 0.d0
            enddo
            enddo
          enddo
       else
        B= 0.d0
        X= 0.d0
        D= 0.d0
        ALUG= 0.d0
      endif

      WW= 0.d0

      if (FTflag.eq.1) then
          do jv= 1, NHYP
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
                AL(9*k-8)= 0.d0
                AL(9*k-7)= 0.d0
                AL(9*k-6)= 0.d0
                AL(9*k-5)= 0.d0
                AL(9*k-4)= 0.d0
                AL(9*k-3)= 0.d0
                AL(9*k-2)= 0.d0
                AL(9*k-1)= 0.d0
                AL(9*k  )= 0.d0
              enddo

              jsU= indexU(j-1)+1
              jeU= indexU(j)
              do k= jsU, jeU
                AU(9*k-8)= 0.d0
                AU(9*k-7)= 0.d0
                AU(9*k-6)= 0.d0
                AU(9*k-5)= 0.d0
                AU(9*k-4)= 0.d0
                AU(9*k-3)= 0.d0
                AU(9*k-2)= 0.d0
                AU(9*k-1)= 0.d0
                AU(9*k  )= 0.d0
              enddo
            enddo
            enddo
          enddo
       else
        AL= 0.d0
        AU= 0.d0
      endif

!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
      valA=            POISSON /      (1.d0-POISSON)
      valB= (1.d0-2.d0*POISSON)/(2.d0*(1.d0-POISSON))
      valX= ELAST*(1.d0-POISSON)/((1.d0+POISSON)*(1.d0-2.d0*POISSON))

      valA= valA * valX
      valB= valB * valX

      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2
        QP1= 1.d0 + POS(ip)
        QM1= 1.d0 - POS(ip)
        EP1= 1.d0 + POS(jp)
        EM1= 1.d0 - POS(jp)
        TP1= 1.d0 + POS(kp)
        TM1= 1.d0 - POS(kp)
        SHAPE(ip,jp,kp,1)= O8th * QM1 * EM1 * TM1
        SHAPE(ip,jp,kp,2)= O8th * QP1 * EM1 * TM1
        SHAPE(ip,jp,kp,3)= O8th * QP1 * EP1 * TM1
        SHAPE(ip,jp,kp,4)= O8th * QM1 * EP1 * TM1
        SHAPE(ip,jp,kp,5)= O8th * QM1 * EM1 * TP1
        SHAPE(ip,jp,kp,6)= O8th * QP1 * EM1 * TP1
        SHAPE(ip,jp,kp,7)= O8th * QP1 * EP1 * TP1
        SHAPE(ip,jp,kp,8)= O8th * QM1 * EP1 * TP1
        PNQ(jp,kp,1)= - O8th * EM1 * TM1
        PNQ(jp,kp,2)= + O8th * EM1 * TM1
        PNQ(jp,kp,3)= + O8th * EP1 * TM1
        PNQ(jp,kp,4)= - O8th * EP1 * TM1
        PNQ(jp,kp,5)= - O8th * EM1 * TP1
        PNQ(jp,kp,6)= + O8th * EM1 * TP1
        PNQ(jp,kp,7)= + O8th * EP1 * TP1
        PNQ(jp,kp,8)= - O8th * EP1 * TP1
        PNE(ip,kp,1)= - O8th * QM1 * TM1
        PNE(ip,kp,2)= - O8th * QP1 * TM1
        PNE(ip,kp,3)= + O8th * QP1 * TM1
        PNE(ip,kp,4)= + O8th * QM1 * TM1
        PNE(ip,kp,5)= - O8th * QM1 * TP1
        PNE(ip,kp,6)= - O8th * QP1 * TP1
        PNE(ip,kp,7)= + O8th * QP1 * TP1
        PNE(ip,kp,8)= + O8th * QM1 * TP1
        PNT(ip,jp,1)= - O8th * QM1 * EM1
        PNT(ip,jp,2)= - O8th * QP1 * EM1
        PNT(ip,jp,3)= - O8th * QP1 * EP1
        PNT(ip,jp,4)= - O8th * QM1 * EP1
        PNT(ip,jp,5)= + O8th * QM1 * EM1
        PNT(ip,jp,6)= + O8th * QP1 * EM1
        PNT(ip,jp,7)= + O8th * QP1 * EP1
        PNT(ip,jp,8)= + O8th * QM1 * EP1
      enddo
      enddo
      enddo

      call MPI_Barrier (SOLVER_COMM, ierr)
      S2_time= MPI_WTIME()
      do icol= 1, ELMCOLORtot
!$omp parallel do private (icel0,icel,in1,in2,in3,in4,in5,in6,in7,in8)  &
!$omp&            private (in10,in20,in30,in40,in50,in60,in70,in80)     &
!$omp&            private (nodLOCAL,ie,je,ip,jp,kk,iiS,iiE,iDlu,k)      &
!$omp&            private (PNXi,PNYi,PNZi,PNXj,PNYj,PNZj,a11,a12)       &
!$omp&            private (a13,a21,a22,a23,a31,a32,a33,ipn,jpn,kpn,coef)&
!$omp&            private (DETJ,PNX,PNY,PNZ,X1,X2,X3,X4,X5,X6,X7,X8)    &
!$omp&            private (Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                     &
!$omp&            private (Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)
      do icel0= ELMCOLORindex(icol-1)+1, ELMCOLORindex(icol)
        icel= ELMCOLORitem(icel0)

!        write (*,*) icol,icel0,icel

        in10= ICELNOD(icel,1)
        in20= ICELNOD(icel,2)
        in30= ICELNOD(icel,3)
        in40= ICELNOD(icel,4)
        in50= ICELNOD(icel,5)
        in60= ICELNOD(icel,6)
        in70= ICELNOD(icel,7)
        in80= ICELNOD(icel,8)

        in1= OLDtoNEW(in10)
        in2= OLDtoNEW(in20)
        in3= OLDtoNEW(in30)
        in4= OLDtoNEW(in40)
        in5= OLDtoNEW(in50)
        in6= OLDtoNEW(in60)
        in7= OLDtoNEW(in70)
        in8= OLDtoNEW(in80)

!C
!C
!C== JACOBIAN & INVERSE JACOBIAN
        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8

        X1= 0.d0
        Y1= 0.d0
        Z1= 0.d0

        X2= DX
        Y2= 0.d0
        Z2= 0.d0

        X3= DX
        Y3= DY
        Z3= 0.d0

        X4= 0.d0
        Y4= DY
        Z4= 0.d0

        X5= 0.d0
        Y5= 0.d0
        Z5= DZ

        X6= DX
        Y6= 0.d0
        Z6= DZ

        X7= DX
        Y7= DY
        Z7= DZ

        X8= 0.d0
        Y8= DY
        Z8= DZ

        call JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,                &
     &               X1, X2, X3, X4, X5, X6, X7, X8,                    &
     &               Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,                    &
     &               Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )

!C
!C== CONSTRUCT the GLOBAL MATRIX
        do ie= 1, 8
          ip = nodLOCAL(ie)
        if (ip.le.N) then
        do je= 1, 8
          jp = nodLOCAL(je)

          IDlu= 0
          if (ip.eq.jp) IDlu= 0

          kk= 0
            iiS= indexU(ip-1) + 1
            iiE= indexU(ip  )
            do k= iiS, iiE
              if ( itemU(k).eq.jp ) then
                kk  = k
                IDlu= 1
                exit
              endif
            enddo

          if (kk.eq.0) then
            iiS= indexL(ip-1) + 1
            iiE= indexL(ip  )

            do k= iiS, iiE
              if ( itemL(k).eq.jp) then
                kk= k
                IDlu= -1
              endif
            enddo
          endif

          PNXi= 0.d0
          PNYi= 0.d0
          PNZi= 0.d0
          PNXj= 0.d0
          PNYj= 0.d0
          PNZj= 0.d0
          
          a11= 0.d0
          a21= 0.d0
          a31= 0.d0
          a12= 0.d0
          a22= 0.d0
          a32= 0.d0
          a13= 0.d0
          a23= 0.d0
          a33= 0.d0
          do kpn= 1, 2
          do jpn= 1, 2
          do ipn= 1, 2
            coef= dabs(DETJ(ipn,jpn,kpn))*WEI(ipn)*WEI(jpn)*WEI(kpn)

            PNXi= PNX(ipn,jpn,kpn,ie)
            PNYi= PNY(ipn,jpn,kpn,ie)
            PNZi= PNZ(ipn,jpn,kpn,ie)

            PNXj= PNX(ipn,jpn,kpn,je)
            PNYj= PNY(ipn,jpn,kpn,je)
            PNZj= PNZ(ipn,jpn,kpn,je)

            a11= a11 + (valX*PNXi*PNXj+valB*(PNYi*PNYj+PNZi*PNZj))*coef
            a22= a22 + (valX*PNYi*PNYj+valB*(PNZi*PNZj+PNXi*PNXj))*coef
            a33= a33 + (valX*PNZi*PNZj+valB*(PNXi*PNXj+PNYi*PNYj))*coef

            a12= a12 + (valA*PNXi*PNYj + valB*PNXj*PNYi)*coef
            a13= a13 + (valA*PNXi*PNZj + valB*PNXj*PNZi)*coef
            a21= a21 + (valA*PNYi*PNXj + valB*PNYj*PNXi)*coef
            a23= a23 + (valA*PNYi*PNZj + valB*PNYj*PNZi)*coef
            a31= a31 + (valA*PNZi*PNXj + valB*PNZj*PNXi)*coef
            a32= a32 + (valA*PNZi*PNYj + valB*PNZj*PNYi)*coef
          enddo
          enddo
          enddo

!          write (*,'(10i8)') icol,icel0,icel,ip,jp,IDlu,kk
          if (IDlu.eq.1) then
            AU(9*kk-8)= AU(9*kk-8) + a11
            AU(9*kk-7)= AU(9*kk-7) + a12
            AU(9*kk-6)= AU(9*kk-6) + a13
            AU(9*kk-5)= AU(9*kk-5) + a21
            AU(9*kk-4)= AU(9*kk-4) + a22
            AU(9*kk-3)= AU(9*kk-3) + a23
            AU(9*kk-2)= AU(9*kk-2) + a31
            AU(9*kk-1)= AU(9*kk-1) + a32
            AU(9*kk  )= AU(9*kk  ) + a33
          endif

          if (IDlu.eq.-1) then
            AL(9*kk-8)= AL(9*kk-8) + a11
            AL(9*kk-7)= AL(9*kk-7) + a12
            AL(9*kk-6)= AL(9*kk-6) + a13
            AL(9*kk-5)= AL(9*kk-5) + a21
            AL(9*kk-4)= AL(9*kk-4) + a22
            AL(9*kk-3)= AL(9*kk-3) + a23
            AL(9*kk-2)= AL(9*kk-2) + a31
            AL(9*kk-1)= AL(9*kk-1) + a32
            AL(9*kk  )= AL(9*kk  ) + a33
          endif

          if (IDlu.eq.0) then
            D(9*ip-8)= D(9*ip-8) + a11
            D(9*ip-7)= D(9*ip-7) + a12
            D(9*ip-6)= D(9*ip-6) + a13
            D(9*ip-5)= D(9*ip-5) + a21
            D(9*ip-4)= D(9*ip-4) + a22
            D(9*ip-3)= D(9*ip-3) + a23
            D(9*ip-2)= D(9*ip-2) + a31
            D(9*ip-1)= D(9*ip-1) + a32
            D(9*ip  )= D(9*ip  ) + a33
          endif

        enddo
        endif
        enddo
      enddo
      enddo
      S3_time= MPI_WTIME()


      if (my_rank.eq.0) write (*,'(1pe16.6,a)') X1_time-S1_time,
     &                         '   coloring of elements'
      if (my_rank.eq.0) write (*,'(1pe16.6,a)') S3_time-S2_time,
     &                         '   matrix assembling'
      
      return
      end
