!********************************************
! N*Nの行列行列積 (C = A * B) にBLASを用いる例
!********************************************

program main
  use omp_lib
  use cublas
  implicit none
  integer :: nt = 100
  integer :: n = 4192
  real(KIND=8),allocatable,dimension(:,:) :: A, B, C, ANS
  real(KIND=8) :: alpha, beta
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: err,dummy(1)
  integer :: i,j,times,ierr
  
  alpha=1.0d0
  beta =0.0d0
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
  ! 一番手っ取り早いのは、OpenACCじゃなくてもOpenACCのdata指示文を使う方法。
  ! Unified Memory使用時は不要なはずだが、一応書いておくとよい
  !******************************
  !$acc data copyin(A,B) copy(C)
  
  !******************************
  ! 実行時間計測時の注意点。
  ! BLASの一回目はやたら時間がかかるようなので除外
  !******************************
  call cublasDgemm('n','n', n, n, n, alpha, A, n, B, n, beta, C, n)
  t2 = omp_get_wtime()
  do concurrent (j=1:n)
     do concurrent (i=1:n)
        C(i,j) = 0.0d0
     end do
  end do

  !******************************
  ! 実行時間計測時の注意点。
  ! 計測する上で必要な同期関数がdo concurretには無いので、
  ! OpenACCのwait指示文を使用
  !******************************   
  t3 = omp_get_wtime()
  call cublasDgemm('n','n', n, n, n, alpha, A, n, B, n, beta, C, n)
  !$acc wait
  t4 = omp_get_wtime()
  !$acc end data

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
  print *, "First run [s]: ", (t2-t1)-(t4-t3)
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

