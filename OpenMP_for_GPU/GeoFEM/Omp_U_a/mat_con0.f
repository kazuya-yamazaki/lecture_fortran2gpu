!C
!C***
!C*** MAT_CON0
!C***
!C
      subroutine MAT_CON0 (ITERkk)
      use hpcmw_fem_mesh
      use hpcmw_fem_util
      use hpcmw_solver_matrix
      implicit REAL*8 (A-H,O-Z)

      N2= 256
      NU= 26
      NL= 26

      if (allocated(INL))  deallocate (INL)
      if (allocated(IAL))  deallocate (IAL)
      if (allocated(INU))  deallocate (INU)
      if (allocated(IAU))  deallocate (IAU)

      allocate (INL(NP), IAL(NP,NL))
      allocate (INU(NP), IAU(NP,NU))

      INL= 0
      IAL= 0
      INU= 0
      IAU= 0

      do icel= 1, ICELTOT
        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

        call FIND_TS_NODE (in1,in2)
        call FIND_TS_NODE (in1,in3)
        call FIND_TS_NODE (in1,in4)
        call FIND_TS_NODE (in1,in5)
        call FIND_TS_NODE (in1,in6)
        call FIND_TS_NODE (in1,in7)
        call FIND_TS_NODE (in1,in8)

        call FIND_TS_NODE (in2,in1)
        call FIND_TS_NODE (in2,in3)
        call FIND_TS_NODE (in2,in4)
        call FIND_TS_NODE (in2,in5)
        call FIND_TS_NODE (in2,in6)
        call FIND_TS_NODE (in2,in7)
        call FIND_TS_NODE (in2,in8)

        call FIND_TS_NODE (in3,in1)
        call FIND_TS_NODE (in3,in2)
        call FIND_TS_NODE (in3,in4)
        call FIND_TS_NODE (in3,in5)
        call FIND_TS_NODE (in3,in6)
        call FIND_TS_NODE (in3,in7)
        call FIND_TS_NODE (in3,in8)

        call FIND_TS_NODE (in4,in1)
        call FIND_TS_NODE (in4,in2)
        call FIND_TS_NODE (in4,in3)
        call FIND_TS_NODE (in4,in5)
        call FIND_TS_NODE (in4,in6)
        call FIND_TS_NODE (in4,in7)
        call FIND_TS_NODE (in4,in8)

        call FIND_TS_NODE (in5,in1)
        call FIND_TS_NODE (in5,in2)
        call FIND_TS_NODE (in5,in3)
        call FIND_TS_NODE (in5,in4)
        call FIND_TS_NODE (in5,in6)
        call FIND_TS_NODE (in5,in7)
        call FIND_TS_NODE (in5,in8)

        call FIND_TS_NODE (in6,in1)
        call FIND_TS_NODE (in6,in2)
        call FIND_TS_NODE (in6,in3)
        call FIND_TS_NODE (in6,in4)
        call FIND_TS_NODE (in6,in5)
        call FIND_TS_NODE (in6,in7)
        call FIND_TS_NODE (in6,in8)

        call FIND_TS_NODE (in7,in1)
        call FIND_TS_NODE (in7,in2)
        call FIND_TS_NODE (in7,in3)
        call FIND_TS_NODE (in7,in4)
        call FIND_TS_NODE (in7,in5)
        call FIND_TS_NODE (in7,in6)
        call FIND_TS_NODE (in7,in8)

        call FIND_TS_NODE (in8,in1)
        call FIND_TS_NODE (in8,in2)
        call FIND_TS_NODE (in8,in3)
        call FIND_TS_NODE (in8,in4)
        call FIND_TS_NODE (in8,in5)
        call FIND_TS_NODE (in8,in6)
        call FIND_TS_NODE (in8,in7)
      enddo

      do in= 1, NP
        NN= INL(in)
        do k= 1, NN
          NCOL1(k)= IAL(in,k)
        enddo
          call mSORT (NCOL1, NCOL2, NN)
        do k= NN, 1, -1
          IAL(in,NN-k+1)= NCOL1(NCOL2(k))
        enddo        

        NN= INU(in)
        do k= 1, NN
          NCOL1(k)= IAU(in,k)
        enddo
          call mSORT (NCOL1, NCOL2, NN)
        do k= NN, 1, -1
          IAU(in,NN-k+1)= NCOL1(NCOL2(k))
        enddo        
      enddo

      contains
!C
!C***
!C*** FIND_TS_NODE
!C***
!C
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
      end subroutine MAT_CON0
