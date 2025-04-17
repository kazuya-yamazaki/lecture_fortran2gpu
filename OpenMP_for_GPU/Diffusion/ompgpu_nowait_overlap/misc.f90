module misc
  implicit none
  real(KIND=8) :: t_s
  real(KIND=8) :: dummy(1)

contains

  subroutine swap(f, fn)
    real(KIND=8),pointer,dimension(:,:,:),intent(inout) :: f,fn
    real(KIND=8),pointer,dimension(:,:,:) :: ftmp

    ftmp => f
    f => fn
    fn => ftmp
  end subroutine swap

  subroutine start_timer()
    real(KIND=8) :: omp_get_wtime
    t_s = omp_get_wtime()
  end subroutine start_timer

  real(KIND=8) function get_elapsed_time()
    real(KIND=8) :: omp_get_wtime
    get_elapsed_time = omp_get_wtime() - t_s
  end function get_elapsed_time
  
  subroutine sync_omp_target()
    integer :: i
    
    !$omp target loop
    do i = 1, 1
      dummy(i) = 0
    enddo
  end subroutine sync_omp_target

end module misc
