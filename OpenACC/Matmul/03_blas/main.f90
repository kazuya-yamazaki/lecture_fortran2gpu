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
  integer :: times
  
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

  !******************************
  ! 実行時間計測時の注意点。
  ! BLASの一回目はやたら時間がかかるようなので除外
  !******************************
  call cublasDgemm('n','n', n, n, n, alpha, A, n, B, n, beta, C, n)
  !$acc kernels
  C(:,:) = 0
  !$acc end kernels
  t3 = omp_get_wtime()
  call cublasDgemm('n','n', n, n, n, alpha, A, n, B, n, beta, C, n)
  !$acc wait
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
    do j = 1, N
       do i = 1, N
          err = err + (a(i,j) - b(i,j)) * (a(i,j) - b(i,j))
       end do
    end do
    accuracy = dble(sqrt(err/(N*N)))
  end function accuracy
  
end program main

   
