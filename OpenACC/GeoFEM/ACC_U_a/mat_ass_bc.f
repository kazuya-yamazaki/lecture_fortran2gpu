!C
!C***
!C*** MAT_ASS_BC
!C***
!C
      subroutine MAT_ASS_BC
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)

      if (allocated(IWKX)) deallocate (IWKX)
      allocate (IWKX(NP,2))
      IWKX= 0

!C
!C== Z=Zmax

!$acc data create(IWKX)                                                                                                                           
!$acc& copyin(BC_STACK, BC_NOD)                                                                                                                   
!$acc& present(OLDtoNEW,ELMCOLORindex,ELMCOLORitem,ICELNOD)                                                                                       
!$acc& present(B,D,AL,AU,indexL,indexU,itemL,itemU)         

!$acc kernels
!$acc loop independent
      do in= 1, NP
        IWKX(in,1)= 0
      enddo
!$acc end kernels
      
      ib0= 4
!$acc kernels
!$acc loop independent
      do ib= BC_STACK(ib0-1)+1, BC_STACK(ib0)
        in= OLDtoNEW(BC_NOD(ib))
        IWKX(in,1)= 1
      enddo
!$acc end kernels
      
      STRESS= ELAST * STRAIN
      VAL   = 0.25d0*DX*DY * STRESS
      do icol= 1, ELMCOLORtot

!$acc kernels         
!$acc loop independent
      do icel0= ELMCOLORindex(icol-1)+1, ELMCOLORindex(icol)
        icel= ELMCOLORitem(icel0)

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

        iq1= IWKX(in1,1)
        iq2= IWKX(in2,1)
        iq3= IWKX(in3,1)
        iq4= IWKX(in4,1)
        iq5= IWKX(in5,1)
        iq6= IWKX(in6,1)
        iq7= IWKX(in7,1)
        iq8= IWKX(in8,1)

        if ((iq1+iq2+iq3+iq4).eq.4) then
          B(3*in10)=  B(3*in10) + VAL
          B(3*in20)=  B(3*in20) + VAL
          B(3*in30)=  B(3*in30) + VAL
          B(3*in40)=  B(3*in40) + VAL
        endif
        if ((iq5+iq6+iq7+iq8).eq.4) then
          B(3*in50)=  B(3*in50) + VAL
          B(3*in60)=  B(3*in60) + VAL
          B(3*in70)=  B(3*in70) + VAL
          B(3*in80)=  B(3*in80) + VAL
        endif
        if ((iq1+iq2+iq5+iq6).eq.4) then
          B(3*in10)=  B(3*in10) + VAL
          B(3*in20)=  B(3*in20) + VAL
          B(3*in50)=  B(3*in50) + VAL
          B(3*in60)=  B(3*in60) + VAL
        endif
        if ((iq2+iq3+iq6+iq7).eq.4) then
          B(3*in20)=  B(3*in20) + VAL
          B(3*in30)=  B(3*in30) + VAL
          B(3*in60)=  B(3*in60) + VAL
          B(3*in70)=  B(3*in70) + VAL
        endif
        if ((iq3+iq4+iq7+iq8).eq.4) then
          B(3*in30)=  B(3*in30) + VAL
          B(3*in40)=  B(3*in40) + VAL
          B(3*in70)=  B(3*in70) + VAL
          B(3*in80)=  B(3*in80) + VAL
        endif
        if ((iq4+iq1+iq8+iq5).eq.4) then
          B(3*in40)=  B(3*in40) + VAL
          B(3*in10)=  B(3*in10) + VAL
          B(3*in80)=  B(3*in80) + VAL
          B(3*in50)=  B(3*in50) + VAL
        endif
      enddo
!$acc end kernels
      enddo

!C
!C== Z=Zmin
!$acc kernels
!$acc loop independent
      do in= 1, NP
        IWKX(in,1)= 0
      enddo
!$acc end kernels
      
      ib0= 3
!$acc kernels
!$acc loop independent
      do ib= BC_STACK(ib0-1)+1, BC_STACK(ib0)
        in= OLDtoNEW(BC_NOD(ib))
        IWKX(in,1)= 1
      enddo
!$acc end kernels
      
!$acc kernels
!$acc loop independent
      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-6)= 0.d0
          D(9*in-3)= 0.d0
          D(9*in-2)= 0.d0
          D(9*in-1)= 0.d0
          D(9*in  )= 1.d0
          B(3*NEWtoOLD(in))= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-2)= 0.d0
            AL(9*k-1)= 0.d0
            AL(9*k  )= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-2)= 0.d0
            AU(9*k-1)= 0.d0
            AU(9*k  )= 0.d0
          enddo
        endif
      enddo
!$acc end kernels
      
!$acc kernels
!$acc loop independent
      do in= 1, N
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-6)= 0.d0
            AL(9*k-3)= 0.d0
            AL(9*k  )= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-6)= 0.d0
            AU(9*k-3)= 0.d0
            AU(9*k  )= 0.d0
          endif
        enddo
      enddo
!$acc end kernels
      
!C
!C== X= Xmin
!$acc kernels
!$acc loop independent
      do in= 1, NP
        IWKX(in,1)= 0
      enddo
!$acc end kernels
      
      ib0= 1
!$acc kernels
!$acc loop independent
      do ib= BC_STACK(ib0-1)+1, BC_STACK(ib0)
        in= OLDtoNEW(BC_NOD(ib))
        IWKX(in,1)= 1
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-8)= 1.d0
          D(9*in-7)= 0.d0
          D(9*in-6)= 0.d0
          D(9*in-5)= 0.d0
          D(9*in-2)= 0.d0
          B(3*NEWtoOLD(in)-2)= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-8)= 0.d0
            AL(9*k-7)= 0.d0
            AL(9*k-6)= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-8)= 0.d0
            AU(9*k-7)= 0.d0
            AU(9*k-6)= 0.d0
          enddo
        endif
      enddo
!$acc end kernels
      
!$acc kernels
!$acc loop independent
      do in= 1, N
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-8)= 0.d0
            AL(9*k-5)= 0.d0
            AL(9*k-2)= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-8)= 0.d0
            AU(9*k-5)= 0.d0
            AU(9*k-2)= 0.d0
          endif
        enddo
      enddo
!$acc end kernels
      
!C
!C== Y= Ymin
!$acc kernels
!$acc loop independent
      do in= 1, NP
        IWKX(in,1)= 0
      enddo
!$acc end kernels
      
      ib0= 2
!$acc kernels
!$acc loop independent
      do ib= BC_STACK(ib0-1)+1, BC_STACK(ib0)
        in= OLDtoNEW(BC_NOD(ib))
        IWKX(in,1)= 1
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-7)= 0.d0
          D(9*in-4)= 1.d0
          D(9*in-1)= 0.d0
          D(9*in-5)= 0.d0
          D(9*in-3)= 0.d0
          B(3*NEWtoOLD(in)-1)= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-5)= 0.d0
            AL(9*k-4)= 0.d0
            AL(9*k-3)= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-5)= 0.d0
            AU(9*k-4)= 0.d0
            AU(9*k-3)= 0.d0
          enddo
        endif
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do in= 1, N
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-7)= 0.d0
            AL(9*k-4)= 0.d0
            AL(9*k-1)= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-7)= 0.d0
            AU(9*k-4)= 0.d0
            AU(9*k-1)= 0.d0
          endif
        enddo
      enddo
!$acc end kernels

!$acc end data
      deallocate (IWKX)

!$acc update host(AL)                                                                 
!$acc update host(AU)                                                                 
!$acc update host(D)                                                                  
!$acc update host(B)
      
      return
      end
