!********************************************
! atomic命令を利用してhistgramを並列化する例
!********************************************

module test
  implicit none
contains
  
  subroutine histgram(x, y)
    integer,dimension(:),intent(in) :: x
    integer,dimension(:),intent(out) :: y 
    integer :: i,j

    !$omp parallel do private(j)
    do i = 1, size(x)
       j = x(i)
       !$omp atomic 
       y(j) = y(j) + 1
       !$omp end atomic
    end do
    !$omp end parallel do 
    
  end subroutine histgram
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 10000000
  integer :: ny = 100
  integer,allocatable,dimension(:) :: y,x
  real(KIND=8) :: a
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: i,times
  
  allocate(x(n),y(ny))

  y(:) = 0
  do i = 1, n
     call random_number(a)
     x(i) = int(a*100)+1
  end do
  
  t1 = omp_get_wtime()
  do times = 1, nt
     call histgram(x,y)
  end do
  t2 = omp_get_wtime()

  !******************************
  ! 答えはyの総和を出力。n*ntになるはず。
  !******************************
  print *, "ans :", sum(y)
  print *, "Time / iteration [s]: ", (t2-t1)/nt

  deallocate(x,y)

end program main

