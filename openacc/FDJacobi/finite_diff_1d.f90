! Running: ./finite_diff_1d (maybe with ACC_NOTIFY=1)

module double
  integer, parameter :: dp = kind(1.0d0) ! double precission
end module double

module diff_finitas
  use double
  use openacc
  implicit none

contains

  pure function f(x)
    real(dp), intent(in) :: x
    real(dp) :: f
    f = 2 !-12*x*x
  end function f

  function finite_diff_1d(a, b, alpha, beta, n ) result(u)
    real(dp), intent(in) :: a, b, alpha, beta !! Interval pts & boundary conditions
    integer, intent(in) :: n ! Number of interior points
    real(dp), dimension(0:n+1) :: u
    real(dp), dimension(1:n) :: new_u
    integer :: iter, i
    integer, parameter :: niters = 10000 ! Number of time iterations
    real(dp) :: x_i, hx

    ! Set up data
    hx = (b-a)/(n+2)
    do i = 0, n+1
       u(i) = alpha + (beta-alpha)*i/(n+1)
    enddo

    call acc_init( acc_device_nvidia )
    !acc data copyin(u), create(new_u)
    !$acc data copy(u, new_u)
    do iter = 1,niters
       !$acc kernels loop
       do i = 1, n
          x_i = a + hx*i
          ! print *, x_i, f(x_i)
          new_u(i) = 0.5_dp*( u(i-1) + u(i+1) + hx*hx*2 ) !hx*hx*f(x_i) )
       enddo
       !$acc kernels loop
       do i = 1,n
          u(i) = new_u(i)
       end do
    enddo
    !$acc end data
  end function finite_diff_1d

end module diff_finitas

program test
  use diff_finitas
  implicit none

  integer, parameter :: n=1000000
  real(dp), dimension(0:n+1) :: u
  integer::i
  real t0, t1

  call cpu_time(t0)
  u = finite_diff_1d(-1._dp, 1._dp, 0._dp, 0._dp, n)
  call cpu_time(t1)
  if(.false.) then
     print *, "Result:"
     print *, u
  end if
  print *, "max=", maxval(u)
  print '("Time = ", g8.3," seconds")', t1-t0
end program test
