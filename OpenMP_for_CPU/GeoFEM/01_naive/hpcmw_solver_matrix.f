!C
!C***
!C*** HPCmw_solver_matrix
!C***
!C
      module hpcmw_solver_matrix
      use hpcmw_util
!C
!C-- ELEMENT Coloring
      integer(kind=kint) :: ELMCOLORtot
      integer(kind=kint), dimension(:), allocatable :: ELMCOLORindex
      integer(kind=kint), dimension(:), allocatable :: ELMCOLORitem

!C
!C-- MATRIX SCALARs
      integer(kind=kint), parameter :: NP0 =  300000*16
      integer(kind=kint), parameter :: NPL0= 4250000*16
      integer(kind=kint), parameter :: NPU0= 4250000*16

      integer(kind=kint) :: N, NP, N2, NL, NU, NPL, NPU, FTflag
      integer(kind=kint) :: NCOLORtot, NHYP, npLX1, npUX1
      integer(kind=kint) :: NLmax, NUmax
      integer(kind=kint) :: NCOLORk
!C
!C-- MATRIX arrays
      real(kind=kreal), dimension(:)  , allocatable :: WS, WR
      real(kind=kreal), dimension(:,:), allocatable :: WW

      real(kind=kreal), dimension(:), allocatable :: D, B, X, ALUG
      real(kind=kreal), dimension(:), allocatable :: ALUG_L, ALUG_U
      real(kind=kreal), dimension(:), allocatable ::  AL
      real(kind=kreal), dimension(:), allocatable ::  AU

      integer(kind=kint), dimension(:), allocatable:: OLDtoNEW, NEWtoOLD

      integer(kind=kint), dimension(:),  allocatable :: INL, INU
      integer(kind=kint), dimension(:,:),allocatable :: IAL, IAU

      integer(kind=kint), dimension(:),  allocatable :: INLmc, INUmc
      integer(kind=kint), dimension(:,:),allocatable :: IALmc, IAUmc

      integer(kind=kint), dimension(:,:),allocatable :: IWKX, IW0


      integer(kind=kint), dimension(:), allocatable :: IVECT
      integer(kind=kint), dimension(:), allocatable :: IW   

      integer(kind=kint), dimension(:),  allocatable :: OLDtoNEWmc
      integer(kind=kint), dimension(:),  allocatable :: NEWtoOLDmc
      integer(kind=kint), dimension(:,:), allocatable :: STACKmc
      integer(kind=kint), dimension(:),  allocatable  :: STACKmcG
      integer(kind=kint), dimension(:),  allocatable :: PEon, COLORon

      character(len=80) :: LINE

      integer(kind=kint), dimension(:), allocatable :: NEWtoOLD_L
      integer(kind=kint), dimension(:), allocatable :: NEWtoOLD_U
      integer(kind=kint), dimension(:), allocatable :: OLDtoNEW_L
      integer(kind=kint), dimension(:), allocatable :: OLDtoNEW_U
      integer(kind=kint), dimension(:), allocatable :: LtoU

      integer(kind=kint), dimension(:), allocatable :: indexL
      integer(kind=kint), dimension(:), allocatable :: indexU
      integer(kind=kint), dimension(:), allocatable :: itemL, itemU
      integer(kind=kint), dimension(:), allocatable :: NLmaxHYP
      integer(kind=kint), dimension(:), allocatable :: NUmaxHYP

      integer(kind=kint), dimension(:),  allocatable :: ICHK
      integer(kind=kint), dimension(:),  allocatable :: IVECmc, IVnew

      end module hpcmw_solver_matrix
