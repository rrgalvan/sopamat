program main
  ! Based on first simple test program in:
  !   https://www.pgroup.com/lit/articles/insider/v4n1a1a.htm
  !
  ! Compiling: pgfortran -o openacc_ex1.exe openacc_ex1.f90 -acc -Minfo
  ! Running: ./openacc_ex1.exe (or ACC_NOTIFY=1 ./openacc_ex1.exe)

  integer :: n        ! size of the vector
  real,dimension(:),allocatable :: a, r, e ! vector, results, expected_results
  integer :: i
  character(10) :: arg1
  if( command_argument_count() .gt. 0 )then
     call get_command_argument( 1, arg1 )
     read(arg1,'(i10)') n
  else
     n = 100000
  endif
  if( n .le. 0 ) n = 100000
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

  ! check the results
  do i = 1,n
     if( r(i) .ne. e(i) )then
        print *, i, r(i), e(i)
        stop 'error found'
     endif
  enddo
  print *, n, 'iterations completed'
end program main
