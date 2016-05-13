program main
  ! Based on first simple test program in:
  !   https://www.pgroup.com/lit/articles/insider/v4n1a1a.htm
  !
  ! Compiling: pgfortran -o openacc_ex1 openacc_ex1.f90 -acc -Minfo
  ! Running: ./openacc_ex1 (or ACC_NOTIFY=1 ./openacc_ex1)

  integer :: n ! size of the vectors
  real, dimension(:), allocatable :: a, r, e ! vector, results, expected_results
  integer :: i
  character(16) :: arg1

  ! Read vectors size from command line:
  if( command_argument_count() .gt. 0 )then
     call get_command_argument( 1, arg1 )
     read(arg1,'(i16)') n
     if( n .le. 0 ) stop "Invalid argument"
  else
     n = 100000
  endif

  ! Initialize vectors:
  allocate(a(n), r(n), e(n))
  do i = 1,n
     a(i) = i*2.0
  enddo

  !$acc kernels loop
  do i = 1,n
     r(i) = a(i) * 2.0
  enddo

  ! Non paralelized loop
  do i = 1,n
     e(i) = a(i) * 2.0
  enddo

  ! Check the results:
  do i = 1,n
     if( r(i) .ne. e(i) )then
        print *, i, r(i), e(i)
        stop 'error found'
     endif
  enddo

  print *, n, 'iterations completed'
end program main
