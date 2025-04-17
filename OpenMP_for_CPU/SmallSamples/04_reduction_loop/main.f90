!********************************************
! 総和など、reduction計算を並列化する例
!********************************************

module test
  implicit none
contains
  
  subroutine reduction(x,sum,n)
    real(KIND=8),dimension(:),intent(in) :: x 
    real(KIND=8),intent(inout) :: sum
    integer,intent(in) :: n
    integer :: i, j, k
    
    !$omp parallel do reduction(+:sum)
    do i = 1, n
       sum = sum + x(i)
    end do
    !$omp end parallel do 

  end subroutine reduction
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000000000
  integer :: n3 = 1000
  real(KIND=8),allocatable,dimension(:) :: x
  real(KIND=8) :: sum
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: times
  
  allocate(x(n))

  x(:) = 1.0d0
  sum = 0.0d0
  
  t1 = omp_get_wtime()
  do times = 1, nt
     call reduction(x,sum,n)
  end do
  t2 = omp_get_wtime()
  
  !******************************
  ! 答えはxの総和nt回足し合わせたもの。
  !******************************
  print *, "ans :", sum
  print *, "Time / iteration [s]: ", (t2-t1)/nt

  deallocate(x)

end program main

