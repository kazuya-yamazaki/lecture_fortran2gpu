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
    ! 多重ループの場合、基本的には各ループにomp loop 
    ! を付与すればおよそ事足りる。
    ! ターゲット指示文を用いる場合、スカラ変数はデフォルトでprivateなので、
    ! OpenMP for CPUのようにprivate(i,j)などは必要ない。
    !***********************************************
    
    !$omp target
    !$omp loop 
    do k = 1, n
       !$omp loop
       do j = 1, n
          !$omp loop 
          do i = 1, n
             y(i,j,k) = a * x(i,j,k) + y(i,j,k)
          end do
       end do
    end do
    !$omp end target
#else
    !***********************************************
    ! OMP 4.xの面倒な点は多重ループの並列化で、
    ! どのループをスレッドブロックレベルで並列化するか、
    ! どのループをスレッドレベルで並列化するか、
    ! 明示的に決めなくてはならない。
    ! 基本的には最外ループをteams distribute,
    ! 最内ループをparallel do で並列化することになる。
    ! ちなみにこのプログラムをCPUで動かすこともできるが、
    ! その場合最外のteams distributeが無視され、
    ! 最内のparallel doがCPUのスレッド並列されることになり、
    ! 全く性能は出ない。（なぜこんな仕様に...）
    !***********************************************
    
    !$omp target
    !$omp teams distribute
    do k = 1, n
       do j = 1, n
          !$omp parallel do 
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

