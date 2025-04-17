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
    ! collapseを駆使して、do concurrentを2つ出現させることで、
    ! コンパイラがスレッドブロックレベルの並列化を外側のdo concurretに、
    ! スレッドレベルの並列化を内側のループに割り当ててくれるよう調整する。
    ! 上手くいくかどうかはコンパイラ次第
    !********************************************

    do concurrent (i = 1:n:BLOCKX,j = 1:n:BLOCKY)
       do k = 1, n, BLOCKX
          do concurrent (ii = 0:BLOCKX-1,jj = 0: BLOCKY-1)
             do kk = 0, BLOCKX-1
                C(i+ii,j+jj) = C(i+ii,j+jj) + A(i+ii,k+kk) * B(k+kk,j+jj)
             end do
          end do
       end do
    end do
    
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
  integer :: i,j,times
  
  allocate(A(n,n),B(n,n),C(n,n),ANS(n,n))

  call random_number(A)
  call random_number(B)
  C = 0.0d0
  ANS = 0.0d0

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

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! ページフォルト（扱おうとしたGPU/CPUのデータがCPU/GPUより古い）が発生した時に
  ! 自動で転送される仕組みなので、データ転送の時間を正確に図るのは難しい。
  ! データ転送はdaxpyのカーネル関数内で発生する。
  ! データ転送時間を除外するなら、以下のように最初の一回目を計測から外す。
  !******************************
  call my_matmul(A,B,C,n)
  t2 = omp_get_wtime()
  do concurrent (j=1:n) local(i)
     do concurrent (i=1:n)
        C(i,j) = 0.0d0
     end do
  end do
  t3 = omp_get_wtime()
  call my_matmul(A,B,C,n)
  t4 = omp_get_wtime()

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
  print *, "First run [s]: ", (t2-t1)
  write (*, '(a,f12.5)'), "elapsed time[sec] : ", t4-t3
  write (*, '(a,f12.5)'), "FLOPS[GFlops]     : ", 2.0D0*n*n*(n-1)/1000000000.0D0/(t4-t3)

  deallocate(A,B,C,ANS)

contains
  double precision function accuracy(a, b)
    double precision,dimension(:,:) :: a,b
    double precision :: err
    integer :: i,j
    err = 0.0
    do concurrent (i = 1: N, j = 1: N) reduce(+:err)
       err = err + (a(i,j) - b(i,j)) * (a(i,j) - b(i,j))
    end do
    accuracy = dble(sqrt(err/(N*N)))
  end function accuracy
  
end program main

