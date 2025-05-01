!********************************************
! 三重ループのDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************
#ifndef OMP
#define OMP 5
#endif

module test
  implicit none
contains
  
  subroutine daxpy3(x, y, a, n)
    real(KIND=8),dimension(:,:,:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:,:,:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i, j, k

#if OMP==5
    !***********************************************
    ! ループ間に他の文が挟まっていないシンプルな多重ループの場合、
    ! loop collapse(ループの数)
    ! とすることで、単なる一重ループとみなすことができる。
    ! スレッド割り当ての効率が改善される可能性がある。
    !***********************************************
    
    !$omp target
    !$omp loop collapse(3)
    do k = 1, n
       do j = 1, n
          do i = 1, n
             y(i,j,k) = a * x(i,j,k) + y(i,j,k)
          end do
       end do
    end do
    !$omp end target
#else
    !***********************************************
    ! OpenMP 4.x では、積極的にcollapseを使うべき。
    ! 3重ループ以上では、collapseを使わないと並列性を活かしきれない。
    !***********************************************
    
    !$omp target
    !$omp teams distribute parallel do collapse(3)
    do k = 1, n
       do j = 1, n
          do i = 1, n
             y(i,j,k) = a * x(i,j,k) + y(i,j,k)
          end do
       end do
    end do
    !$omp end target
#endif    
    
  end subroutine daxpy3
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000
  real(KIND=8),allocatable,dimension(:,:,:) :: y,x
  real(KIND=8) :: a
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: i,j,k,times
  
  allocate(y(n,n,n),x(n,n,n))

  y(:,:,:) = 0.0d0
  x(:,:,:) = 1.0d0
  a = 2.0d0
  
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
  ! DAXPY() NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! 上記のようにCPUで初期化したxやyは、まずCPUメモリに置かれ、
  ! カーネル関数を呼んでからマイグレーションが発生する。
  ! マイグレーションの影響を除外するなら、以下のように最初の一回目を計測から外す。
  ! (実際にはいつマイグレーションが起きるかは処理系任せなので、二回目以降にマイグレーションが起きる可能性もある)
  !******************************
  call daxpy3(x,y,a,n)
  t2 = omp_get_wtime()
  !$omp target
  !$omp loop 
  do k = 1, n
     !$omp loop 
     do j = 1, n
        !$omp loop 
        do i = 1, n
           y(i,j,k) = 0.0d0
        end do
     end do
  end do
  !$omp end target
  t3 = omp_get_wtime()
  do times = 1, nt
     call daxpy3(x,y,a,n)
  end do
  t4 = omp_get_wtime()

  !******************************
  ! 答えはyの平均値を出力。2*ntになるはず。
  ! Unified Memoryを使う場合、sum(y)の計算では、CPUがGPUメモリのyを直接参照する。
  ! Managed Memoryでは、sum(y)の計算中に、yの値がページ単位でGPUメモリからCPUメモリに転送され、CPUはCPUメモリ上のyを参照する。
  !******************************
  ans = sum(y(:,:,:))/(n*n*n)
  print *, "ans :", ans
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

