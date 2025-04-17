!C
!C***
!C*** HPCmw_fem_mesh
!C***
!C
      module hpcmw_fem_mesh
      use hpcmw_util
!C
!C-- COMM. TABLE
      integer          :: NEIBPETOT
      integer, dimension(0:26) :: NOD_STACK_IMPORT
      integer, dimension(0:26) :: NOD_STACK_EXPORT
      integer, dimension(300000*4):: NOD_IMPORT, NOD_EXPORT
      integer, dimension(300000*4):: NOD_EXPORT_NEW
      integer, dimension(26)   :: NEIBPE

!C
!C-- PROBLEM PARAMETERs
      real   (kind=kreal), parameter :: DX= 1.d0
      real   (kind=kreal), parameter :: DY= 1.d0
      real   (kind=kreal), parameter :: DZ= 1.d0
      real   (kind=kreal), parameter :: STRAIN = 1.00d0
      real   (kind=kreal), parameter :: ELAST  = 1.00d0
      real   (kind=kreal), parameter :: POISSON= 0.30d0

!C
!C-- CONNECTIVITIES & BOUNDARY nodes
      integer:: ICELTOT
      integer(kind=kint), dimension(:,:),allocatable :: iCELNOD
      integer(kind=kint), dimension(:),  allocatable :: BC_STACK
      integer(kind=kint), dimension(:),  allocatable :: BC_NOD

      end module hpcmw_fem_mesh
