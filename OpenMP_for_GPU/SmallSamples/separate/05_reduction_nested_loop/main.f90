!********************************************
! 総和など、reduction計算を並列化する例
!********************************************
#ifndef OMP
#define OMP 5
#endif

module test
  implicit none
contains
  
  subroutine reduction(y,sum,n)
    real(KIND=8),dimension(:,:,:),intent(in) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),intent(inout) :: sum
    integer,intent(in) :: n
    integer :: i, j, k

#if OMP==5    
    !***********************************************
    ! 多重ループの場合は、reductionが必要なループ全てに付ける。
    ! collapseを使って一重化しても良い。
    !***********************************************
    
    !$omp target
    !$omp loop reduction(+:sum)
    do k = 1, n
       !$omp loop reduction(+:sum)
       do j = 1, n
          !$omp loop reduction(+:sum)
          do i = 1, n
             sum = sum + y(i,j,k)
          end do
       end do
    end do
    !$omp end target

#else
    !***********************************************
    ! OpenMP 4.x で多重ループの場合は、基本はcollapse。
    ! 下のループはcollapseできるのでしたほうがいいが、
    ! collapseできないループを仮定すると、以下のようになる。
    !***********************************************
    
    !$omp target
    !$omp teams distribute reduction(+:sum)
    do k = 1, n
       do j = 1, n
          !$omp parallel do reduction(+:sum)
          do i = 1, n
             sum = sum + y(i,j,k)
          end do
       end do
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
  !$omp target data map(to:y) 
  !******************************
  ! 上のdata指示文のタイミングで、GPU上に配列x, yが確保される。
  ! mapの中の指示子がtofrom, toならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t3 = omp_get_wtime()
  do times = 1, nt
     call reduction(y,sum,n)
  end do
  t4 = omp_get_wtime()
  !$omp end target data
  !******************************
  ! 上のend target dataのタイミングで、GPU上の配列x, yが解放される。
  ! mapの中の指示子がtofrom, fromならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t2 = omp_get_wtime()
  
  !******************************
  ! 答えはxとyの総和nt回足し合わせたもの。
  !******************************
  print *, "ans :", sum
  print *, "CPU-GPU data transfer [s]: ", (t2-t1)-(t4-t3)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(y)

end program main

