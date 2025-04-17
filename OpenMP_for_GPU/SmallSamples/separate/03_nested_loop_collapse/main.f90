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
  integer :: times
  
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
  !$omp target data map(to:x) map(tofrom:y)
  !******************************
  ! 上のdata指示文のタイミングで、GPU上に配列x, yが確保される。
  ! mapの中の指示子がtofrom, toならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t3 = omp_get_wtime()
  do times = 1, nt
     call daxpy3(x,y,a,n)
  end do
  t4 = omp_get_wtime()
  !$omp end target data
  !******************************
  ! 上のend target dataのタイミングで、GPU上の配列x, yが解放される。
  ! mapの中の指示子がtofrom, fromならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t2 = omp_get_wtime()

  
  !******************************
  ! 答えはyの平均値を出力。2*ntになるはず。
  !******************************
  ans = sum(y(:,:,:))/(n*n*n)
  print *, "ans :", ans
  print *, "CPU-GPU data transfer [s]: ", (t2-t1)-(t4-t3)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main
