!********************************************
! 総和など、reduction計算を並列化する例
!********************************************
#ifndef OMP
#define OMP 5
#endif

module test
  implicit none
contains
  
  subroutine reduction(x,sum,n)
    real(KIND=8),dimension(:),intent(in) :: x ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),intent(inout) :: sum
    integer,intent(in) :: n
    integer :: i, j, k
    
#if OMP==5
    !***********************************************
    ! loop の後に、reduction節を付与。OpenMP for CPUとほとんど変わらない。
    !***********************************************

    !$omp target
    !$omp loop reduction(+:sum)
    do i = 1, n
       sum = sum + x(i)
    end do
    !$omp end target

#else
    !***********************************************
    ! teams ~ の後に、reduction節を付与。OpenMP for CPUとほとんど変わらない。
    !***********************************************

    !$omp target
    !$omp teams distribute parallel do reduction(+:sum)
    do i = 1, n
       sum = sum + x(i)
    end do
    !$omp end target

#endif    

    
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
  
  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! Unified Memoryを用いないOpenMPプログラムでは、target指示文を触った時に表面化する。
  ! 次のループはGPUの起動時間を除外するためだけのもので、プログラム上の意味はない。
  !******************************

  t1 = omp_get_wtime()
  !$omp target 
  dummy(1) = 0
  !$omp end target
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  !******************************
  ! NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! 上記のようにCPUで初期化したxやyは、まずCPUメモリに置かれ、
  ! カーネル関数を呼んでからマイグレーションが発生する。
  ! マイグレーションの影響を除外するなら、以下のように最初の一回目を計測から外す。
  ! (実際にはいつマイグレーションが起きるかは処理系任せなので、二回目以降にマイグレーションが起きる可能性もある)
  !******************************
  call reduction(x,sum,n)
  t2 = omp_get_wtime()
  sum = 0.0d0
  t3 = omp_get_wtime()
  do times = 1, nt
     call reduction(x,sum,n)
  end do
  t4 = omp_get_wtime()
  
  !******************************
  ! 答えはxの総和nt回足し合わせたもの。
  !******************************
  print *, "ans :", sum
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x)

end program main

