!********************************************
! 単純な一重ループとしてDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************

module test
  implicit none
contains

  real(KIND=8) function madd(a,x,y)
    real(KIND=8),intent(in) :: a,x,y
    !$acc routine seq

    !***********************************************
    ! この関数がGPUで実行されることをroutine指示文でコンパイラに教える。
    ! そうしないと、GPU向けのバイナリが生成されない。
    ! gang, (worker), vector, seq を選択できるが、ループの最内側で呼ぶ場合はseq。
    ! gang or vectorを指定された関数はスレッドブロック or スレッドのグループで実行され、
    ! 関数内でさらにloop指示文による並列化も可能（kernelsは付けない）。
    !***********************************************

    madd = a * x + y
  end function madd
  
  subroutine daxpy(x, y, a, n)
    real(KIND=8),dimension(:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i
    !$acc routine(madd) 
    
    !***********************************************
    ! 並列化するループ内から関数呼び出しを行う場合、
    ! routine 指示文が必要となる。
    ! 呼び出し側でも、上のようにどの関数がGPU実行可能なのかを示す必要がある。
    !***********************************************
    
    !$acc kernels
    !$acc loop independent
    do i = 1, n
       y(i) = madd(a,x(i),y(i))
    end do
    !$acc end kernels
    
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
  integer :: times
  
  allocate(y(n),x(n))

  y(:) = 0.0d0
  x(:) = 1.0d0
  a = 2.0d0
  
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
  call daxpy(x,y,a,n)
  t2 = omp_get_wtime()
  !$acc kernels
  y(:) = 0.0d0
  !$acc end kernels
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

