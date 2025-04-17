!********************************************
! 単純な一重ループとしてDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************

module test
  implicit none
contains
  
  subroutine daxpy(x, y, a, n)
    real(KIND=8),dimension(:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i
    
    !***********************************************
    ! 基本的にはdo を do concurrentに置き換えるだけ
    !***********************************************
    
    do concurrent (i = 1: n)
       y(i) = a * x(i) + y(i)
    end do
    
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
  real(KIND=8) :: t1,t2,t3,t4,t5,t6
  real(KIND=8) :: ans, dummy(1)
  integer :: i,times
  
  allocate(y(n),x(n))
  
  !******************************
  ! Unified Memoryでは、変数の初期化をCPUでやりたい場合でも、
  ! 変数が大きい場合はまずGPUでダミーデータを代入しておくことで、
  ! 変数が最初からGPUメモリに置かれ、マイグレーションのオーバーヘッドを回避できる。
  ! Managed Memoryを用いる場合、GPUでのダミーデータ代入は初期化が遅くなるだけで無意味。
  !******************************
#if defined(__NVCOMPILER_GPU_UNIFIED_MEM) && !defined(__NVCOMPILER_GPU_MANAGED_MEM)
  do concurrent (i=1:n)
     y(i) = 0.0d0
     x(i) = 0.0d0
  end do
#endif

  y(:) = 0.0d0
  x(:) = 1.0d0
  a = 2.0d0

  !******************************
  ! DAXPY() NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! 上記のようなダミーデータ代入を行わなかった場合は、
  ! xやyはまずCPUメモリに置かれ、daxpyのカーネル関数を呼んでからマイグレーションが発生する。
  ! マイグレーションの影響を除外するなら、以下のように最初の一回目を計測から外す。
  ! (実際にはいつマイグレーションが起きるかは処理系任せなので、二回目以降にマイグレーションが起きる可能性もある)
  !******************************
  call daxpy(x,y,a,n)
  t2 = omp_get_wtime()
  do concurrent (i=1:n)
     y(i) = 0.0d0
  end do
  t3 = omp_get_wtime()
  do times = 1, nt
     call daxpy(x,y,a,n)
  end do
  t4 = omp_get_wtime()

  
  !******************************
  ! 答えはyの平均値を出力。2*ntになるはず。
  ! Unified Memoryを使う場合、sum(y)の計算では、CPUがGPUメモリのyを直接参照する。
  ! Managed Memoryでは、sum(y)の計算中に、yの値がページ単位でGPUメモリからCPUメモリに転送され、CPUはCPUメモリ上のyを参照する。
  !******************************
  ans = sum(y)/n
  print *, "ans :", ans
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

