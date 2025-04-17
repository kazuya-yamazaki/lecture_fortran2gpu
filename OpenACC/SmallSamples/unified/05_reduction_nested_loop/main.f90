!********************************************
! 総和など、reduction計算を並列化する例
!********************************************

module test
  implicit none
contains
  
  subroutine reduction(y,sum,n)
    real(KIND=8),dimension(:,:,:),intent(in) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),intent(inout) :: sum
    integer,intent(in) :: n
    integer :: i, j, k
    
    !***********************************************
    ! 多重ループの場合は、reductionが必要なループ全てに付ける。
    ! collapseを使って一重化しても良い。
    !***********************************************
    
    !$acc kernels
    !$acc loop independent reduction(+:sum)
    do k = 1, n
       !$acc loop independent reduction(+:sum)
       do j = 1, n
          !$acc loop independent reduction(+:sum)
          do i = 1, n
             sum = sum + y(i,j,k)
          end do
       end do
    end do
    !$acc end kernels
    
  end subroutine reduction
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000
  real(KIND=8),allocatable,dimension(:,:,:) :: y
  real(KIND=8) :: sum
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: times
  
  allocate(y(n,n,n))

  y(:,:,:) = 1.0d0
  sum = 0.0d0
  
  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! Unified Memoryを用いる場合、プログラムの開始時にすでに起動処理がかかるので、
  ! プログラム中の最初の指示文では表面化しない。
  ! 下のkernelsにかかる時間は通常のOpenACCより小さい
  !******************************

  t1 = omp_get_wtime()
  !$acc kernels
  dummy(:) = 0
  !$acc end kernels
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  !******************************
  ! NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! 上記のようにCPUで初期化したxやyは、まずCPUメモリに置かれ、
  ! daxpyのカーネル関数を呼んでからマイグレーションが発生する。
  ! マイグレーションの影響を除外するなら、以下のように最初の一回目を計測から外す。
  ! (実際にはいつマイグレーションが起きるかは処理系任せなので、二回目以降にマイグレーションが起きる可能性もある)
  !******************************
  t1 = omp_get_wtime()
  call reduction(y,sum,n)
  t2 = omp_get_wtime()
  sum = 0.0d0
  t3 = omp_get_wtime()
  do times = 1, nt
     call reduction(y,sum,n)
  end do
  t4 = omp_get_wtime()
  
  !******************************
  ! 答えはxとyの総和nt回足し合わせたもの。
  !******************************
  print *, "ans :", sum
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(y)

end program main

