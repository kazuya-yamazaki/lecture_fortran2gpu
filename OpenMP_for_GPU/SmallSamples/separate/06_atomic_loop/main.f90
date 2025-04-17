!********************************************
! atomic命令を利用してhistgramを並列化する例
!********************************************
#ifndef OMP
#define OMP 5
#endif

module test
  implicit none
contains
  
  subroutine histgram(x, y)
    integer,dimension(:),intent(in) :: x
    integer,dimension(:),intent(out) :: y 
    integer :: i,j

#if OMP==5    
    !***********************************************
    ! atomicの書き方はOpenMP for CPU とほとんど同じ。
    ! NVIDIAのハイエンドGPUでは、atomic addに対して
    ! ハードウェアサポートがあるので、それなりに速い。
    !***********************************************
    
    !$omp target
    !$omp loop 
    do i = 1, size(x)
       j = x(i)
       !$omp atomic 
       y(j) = y(j) + 1
       !$omp end atomic
    end do
    !$omp end target

#else
    !***********************************************
    ! OpenMP 4.xでも書き方は変わらない。
    !***********************************************
    
    !$omp target
    !$omp teams distribute parallel do 
    do i = 1, size(x)
       j = x(i)
       !$omp atomic 
       y(j) = y(j) + 1
       !$omp end atomic
    end do
    !$omp end target

#endif
    
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
  !$omp target data map(to:x) map(tofrom:y)
  !******************************
  ! 上のdata指示文のタイミングで、GPU上に配列x, yが確保される。
  ! mapの中の指示子がtofrom, toならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t3 = omp_get_wtime()
  do times = 1, nt
     call histgram(x,y)
  end do
  t4 = omp_get_wtime()
  !$omp end target data
  !******************************
  ! 上のend target dataのタイミングで、GPU上の配列x, yが解放される。
  ! mapの中の指示子がtofrom, fromならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t2 = omp_get_wtime()

  !******************************
  ! 答えはyの総和を出力。n*ntになるはず。
  !******************************
  print *, "ans :", sum(y)
  print *, "CPU-GPU data transfer [s]: ", (t2-t1)-(t4-t3)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

