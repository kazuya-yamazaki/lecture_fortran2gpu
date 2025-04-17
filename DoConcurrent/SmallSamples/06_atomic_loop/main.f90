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
    
    !***********************************************
    ! do concurrentでatomicはサポートされていない(Fortran 2018）
    ! 以下のコードは誤り。
    !***********************************************
    
    do concurrent (i = 1: size(x))
       j = x(i)
       y(j) = y(j) + 1
    end do
    
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
  
  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! do concurrentはUnified Memoryを用いるので、
  ! プログラムの開始時にすでに起動処理がかかり、
  ! 最初のGPU実行で起動時間が表面化しない。
  !******************************

  t1 = omp_get_wtime()
  do concurrent (i=1:1)
     dummy(1) = 0.0d0
  end do
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  !******************************
  ! NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! 上記のようにCPUで初期化したxやyは、まずCPUメモリに置かれ、
  ! daxpyのカーネル関数を呼んでからマイグレーションが発生する。
  ! マイグレーションの影響を除外するなら、以下のように最初の一回目を計測から外す。
  ! (実際にはいつマイグレーションが起きるかは処理系任せなので、二回目以降にマイグレーションが起きる可能性もある)
  !******************************
  call histgram(x,y)
  t2 = omp_get_wtime()
  do concurrent (i = 1:n)
     y(i) = 0
  end do
  t3 = omp_get_wtime()
  do times = 1, nt
     call histgram(x,y)
  end do
  t4 = omp_get_wtime()

  !******************************
  ! 答えはyの総和を出力。n*ntになるはず。
  !******************************
  print *, "ans :", sum(y)
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

