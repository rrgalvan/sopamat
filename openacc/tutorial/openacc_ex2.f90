program main
  ! Based on second test program in:
  !   https://www.pgroup.com/lit/articles/insider/v4n1a1a.htm
  !
  ! Compiling: pgfortran -o openacc_ex2 openacc_ex2.f90 -acc -Minfo
  ! Running: ./openacc_ex2.exe (or ACC_NOTIFY=1 ./openacc_ex2.exe)

  use openacc
  integer :: n        ! size of the vector
  real(kind(1.d0)),dimension(:),allocatable :: a  ! the vector
  real(kind(1.d0)),dimension(:),allocatable :: r  ! the results
  real(kind(1.d0)),dimension(:),allocatable :: e  ! expected results
  integer :: i
  integer :: c0, c1, c2, c3, cgpu, chost
  character(10) :: arg1
  if( command_argument_count() .gt. 0 )then
     call get_command_argument( 1, arg1 )
     read(arg1,'(i10)') n
  else
     n = 100000
  endif
  if( n .le. 0 ) n = 100000
  allocate(a(n))
  allocate(r(n))
  allocate(e(n))
  do i = 1,n
     a(i) = i*2.0
  enddo
  call acc_init( acc_device_nvidia )
!  do j = 1, 3
     call system_clock( count=c1 )
     !$acc kernels loop
     do i = 1,n
        r(i) = sin(a(i)) ** 2 + cos(a(i)) ** 2
     enddo
     call system_clock( count=c2 )
     cgpu = c2 - c1
     do i = 1,n
        e(i) = sin(a(i)) ** 2 + cos(a(i)) ** 2
     enddo
     call system_clock( count=c3 )
     chost = c3 - c2
     ! check the results
     do i = 1,n
        if( abs(r(i) - e(i)) .gt. 0.000001 )then
           print *, i, r(i), e(i)
        endif
     enddo
     print *, n, ' iterations completed'
     print *, cgpu, ' microseconds on GPU'
     print *, chost, ' microseconds on host'
!  end do

end program main
