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
    f = -12*x*x
  end function f

  function finite_diff_1d(a, b, alpha, beta, npoints ) result(u)
    real(dp), intent(in) :: a, b, alpha, beta !! Interval pts & boundary conditions
    integer, intent(in) :: npoints ! Number of interior points
    real(dp), dimension(0:npoints+1) :: u
    real(dp), dimension(1:npoints) :: new_u
    integer :: iter, i
    integer, parameter :: niters = 5 ! Number of time iterations
    real(dp) :: x_i, hx

    ! Set up data
    hx = (b-a)/(npoints+2)
    do i = 0, npoints+1
       u(i) = alpha + (beta-alpha)*i/(npoints+1)
    enddo

    call acc_init( acc_device_nvidia )
    !$acc data copy(u(:),new_u(:))
    do iter = 1,niters
       !acc kernels loop
       do i = 1, npoints
          x_i = a + hx*i
          ! print *, x_i, f(x_i)
          new_u(i) = 0.5_dp*( u(i-1) + u(i+1) + hx*hx*f(x_i) ) ! 12*x_i*x_i )
       enddo
       print *, iter, new_u
       !acc kernels loop
       do i = 1, npoints
          u(i) = new_u(i)
       enddo
    enddo
    !$acc end data
  end function finite_diff_1d

end module diff_finitas

program test
  use diff_finitas
  implicit none

  integer, parameter :: npoints=10
  real(dp), dimension(0:npoints+1) :: u
  integer::i
  real t0, t1

  call cpu_time(t0)
  u = finite_diff_1d(0._dp, 1._dp, 0._dp, 1._dp, npoints)
  call cpu_time(t1)
  if(.true.) then
     print *, "Result:"
     do i = 0, npoints+1, npoints/20
        write(*, '(g8.3)', advance='no') u(i)
     enddo
     print *
  end if
  print '("Time = ", g8.3," seconds")', t1-t0
end program test
