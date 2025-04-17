!********************************************
! N*Nの行列行列積 (C = A * B) にブロッキングを適用する例
! 行列行列積にはBLASライブラリ等を使うべきだが、
! ブロッキングを適用した複雑なループに対応できるのかどうか
! 調べるためのプログラム。
!********************************************

#define BLOCKX 16
#define BLOCKY 16

module test
  implicit none
contains

  subroutine my_matmul(A, B, C, n)
    double precision, dimension(:,:) :: A, B, C
    integer :: n
    integer :: i, j, k, ii, jj, kk
    
    !********************************************
    ! 各ループをどのレベルで並列化するべきかを見極め、パラメータを指定する必要がある。
    ! OpenACCでは基本的に、gang, vectorを用いてスレッド数などを指定する。
    ! gangはスレッドブロックレベル（CPUでいうコア並列に相当）の並列化、
    ! vectorはスレッドレベル（CPUでいうところのSIMD並列に相当）の並列化を行う。
    ! L1キャッシュはvectorレベルでしか共有されない。
    ! 下で用いているtileは、vectorレベルで行列を平面分割したい時に用いる。
    ! gang, vector, tileなどの設定に関しては、OpenMPからの移行初期段階で考える必要はない。
    !********************************************

    !$acc kernels present(A,B,C)
    !$acc loop collapse(2) independent gang
    do j = 1, n, BLOCKY
       do i = 1, n, BLOCKX
          !$acc loop seq
          do k = 1, n, BLOCKX
             !$acc loop independent tile(BLOCKX,BLOCKY)
             do jj = 0, BLOCKY-1
                do ii = 0, BLOCKX-1
                   do kk = 0, BLOCKX-1
                      C(i+ii,j+jj) = C(i+ii,j+jj) + A(i+ii,k+kk) * B(k+kk,j+jj)
                   end do
                end do
             end do
          end do
       end do
    end do
    !$acc end kernels

  end subroutine my_matmul
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 4192
  real(KIND=8),allocatable,dimension(:,:) :: A, B, C, ANS
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: err,dummy(1)
  integer :: times
  
  allocate(A(n,n),B(n,n),C(n,n),ANS(n,n))

  call random_number(A)
  call random_number(B)
  C = 0.0d0
  ANS = 0.0d0

  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! Unified Memoryを用いないOpenACCプログラムでは、最初にデータ指示文や
  ! kernels指示文を触った時に表面化する。
  ! 次のkernels指示文はGPUの起動時間を除外するためだけのもので、プログラム上の意味はない。
  !******************************

  t1 = omp_get_wtime()
  !$acc kernels
  dummy(:) = 0
  !$acc end kernels
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  t1 = omp_get_wtime()
  !$acc data copyin(A,B) copy(C)
  !******************************
  ! 上のacc dataのタイミングで、GPU上に配列x, yが確保される。
  ! 指示子がcopy, copyinならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t3 = omp_get_wtime()
  call my_matmul(A,B,C,n)
  t4 = omp_get_wtime()
  !$acc end data
  !******************************
  ! 上のacc end dataのタイミングで、GPU上の配列x, yが解放される。
  ! 指示子がcopy, copyoutならば、このタイミングでCPUからGPUへデータの転送が行われる。
  !******************************
  t2 = omp_get_wtime()

  !******************************
  ! 答えは組み込み関数のmatmulと比較
  !******************************
  ANS = matmul(A,B)
  err = accuracy(C, ANS)
  if(err < 0.000001) then
     write (*, '(a)'), "check result...OK"
  else
     write (*, '(a)'), "check result...error!"
  end if
  print *, "CPU-GPU data transfer [s]: ", (t2-t1)-(t4-t3)
  write (*, '(a,f12.5)'), "elapsed time[sec] : ", t4-t3
  write (*, '(a,f12.5)'), "FLOPS[GFlops]     : ", 2.0D0*n*n*(n-1)/1000000000.0D0/(t4-t3)

  deallocate(A,B,C,ANS)

contains
  double precision function accuracy(a, b)
    double precision,dimension(:,:) :: a,b
    double precision :: err
    integer :: i,j
    err = 0.0
    !$acc kernels
    !$acc loop collapse(2) reduction(+:err)
    do j = 1, N
       do i = 1, N
          err = err + (a(i,j) - b(i,j)) * (a(i,j) - b(i,j))
       end do
    end do
    !$acc end kernels
    accuracy = dble(sqrt(err/(N*N)))
  end function accuracy
  
end program main

