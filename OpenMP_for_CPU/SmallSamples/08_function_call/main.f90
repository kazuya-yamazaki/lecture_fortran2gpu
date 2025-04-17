!********************************************
! 単純な一重ループとしてDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************

module test
  implicit none
contains

  real(KIND=8) function madd(a,x,y)
    real(KIND=8),intent(in) :: a,x,y
    !$omp declare simd
    madd = a * x + y
  end function madd
  
  !***********************************************
  ! Intel compilerではコンパイラのインライン展開でなんとかなってしまうことが多いので
  ! あまり認知されていないが、SIMDループ内での関数呼び出しを行う場合、
  ! 上で使っているdeclare simd指示文が必要。
  ! なお、富士通コンパイラではサポートされておらず、SIMD化のためには
  ! 最内ループ内での関数呼び出しを全て手動でインライン展開しなくてはならない。
  !***********************************************
  
  subroutine daxpy(x, y, a, n)
    real(KIND=8),dimension(:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i
    
    !$omp parallel do simd
    do i = 1, n
       y(i) = madd(a,x(i),y(i))
    end do
    !$omp end parallel do 
    
  end subroutine daxpy
  
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000000000
  real(KIND=8),allocatable,dimension(:) :: y,x
  real(KIND=8) :: a
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: times
  
  allocate(y(n),x(n))

  y(:) = 0.0d0
  x(:) = 1.0d0
  a = 2.0d0
  
  t1 = omp_get_wtime()
  do times = 1, nt
     call daxpy(x,y,a,n)
  end do
  t2 = omp_get_wtime()

  
  !******************************
  ! 答えはyの平均値を出力。2*ntになるはず。
  !******************************
  ans = sum(y)/n
  print *, "ans :", ans
  print *, "Time / iteration [s]: ", (t2-t1)/nt

  deallocate(x,y)

end program main

