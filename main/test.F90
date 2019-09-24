program test
  use NdSpline
  implicit none
  integer :: n_dim

  write(*,"(a)") "Enter the test dimenstion 1-4:"
  read(*,*) n_dim

  select case(n_dim)

  case(1)
    call test_1d()
  case(2)
    call test_2d()
  case(3)
    call test_3d()
  case(4)
    !call test_4d()
  case default
  end select

contains

  subroutine test_1d()
    type(spline) :: sp
    real(8) :: xmin = 0.d0, xmax = 8.d0
    integer :: n, k = 4
    real(8), allocatable :: x(:), f(:)
    real(8) :: max_err
    integer :: i

    n = 30
    allocate(x(n), f(n))
    do i = 1, n
      x(i) = xmin + dble(i-1) / dble(n-1) * (xmax - xmin)
    end do

    do i = 1, n
      f(i) = f1(x(i))
    end do

    call sp%init([k],[n],x,f)
    deallocate(x, f)

    n = 100
    allocate(x(n))
    do i = 1, n
      x(i) = xmin + dble(i-1) / dble(n-1) * (xmax - xmin)
    end do
    f = sp%interpolate([n], x)

    max_err = 0.d0
    open(15, file="file_1d.dat")
    write(15,*) "x,  interpolated, true value"
    do i = 1, n
      write(15,"(f8.4, 2es18.6)") x(i), f(i), f1(x(i))
      max_err = max(max_err, abs(f(i) - f1(x(i))))
    end do
    close(15)
    write(*,"(a,es18.6)") "Maximum interpolation error: ", max_err
  end subroutine test_1d

  subroutine test_2d()
    type(spline) :: sp
    real(8) :: xmin = 0.d0, xmax = 8.d0
    real(8) :: ymin = 0.d0, ymax = 8.d0
    integer :: nx, ny, kx=4, ky=4
    real(8), allocatable :: x(:), y(:), f(:,:), f1d(:)
    real(8) :: max_err
    integer :: i, j

    nx = 20
    ny = 20
    allocate(x(nx), y(ny), f(nx,ny), f1d(nx*ny))
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do

    do i = 1, nx
      do j = 1, ny
        f(i,j) = f2(x(i),y(j))
      end do
    end do

    f1d = reshape(f, shape(f1d))
    write(*,*) [x,y]
    call sp%init([kx,ky],[nx,ny],[x,y],f1d)
    deallocate(x, y, f, f1d)

    nx = 100
    ny = 100
    allocate(x(nx),y(ny),f(nx,ny))
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do
    f1d = sp%interpolate([nx, ny], [x,y])
    f = reshape(f1d, shape(f))

    max_err = 0.d0
    open(15, file="file_2d.dat")
    write(15,*) "x,  y, interpolated, true value"
    do i = 1, nx
      do j = 1, ny
        write(15,"(2f8.4, 2es18.6)") x(i), y(j), f(i,j), f2(x(i),y(j))
        max_err = max(max_err, abs(f(i,j) - f2(x(i),y(j))))
      end do
    end do
    close(15)
    write(*,"(a,es18.6)") "Maximum interpolation error: ", max_err
  end subroutine test_2d

  subroutine test_3d()
    type(spline) :: sp
    real(8) :: xmin = 0.d0, xmax = 2.d0
    real(8) :: ymin = 0.d0, ymax = 4.d0
    real(8) :: zmin = 0.d0, zmax = 8.d0
    integer :: nx, ny, nz, kx=4, ky=4, kz = 4
    real(8), allocatable :: x(:), y(:), z(:), f(:,:,:), f1d(:)
    real(8) :: max_err
    integer :: i, j, k

    nx = 10
    ny = 20
    nz = 30
    allocate(x(nx), y(ny), z(nz), f(nx,ny,nz), f1d(nx*ny*nz))
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do
    do i = 1, nz
      z(i) = zmin + dble(i-1) / dble(nz-1) * (zmax - zmin)
    end do

    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          f(i,j,k) = f3(x(i),y(j),z(k))
        end do
      end do
    end do

    f1d = reshape(f, shape(f1d))
    call sp%init([kx,ky,kz],[nx,ny,nz],[x,y,z],f1d)
    deallocate(x, y, z, f, f1d)

    nx = 100
    ny = 100
    nz = 100
    allocate(x(nx),y(ny),z(nz),f(nx,ny,nz))
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do
    do i = 1, nz
      z(i) = zmin + dble(i-1) / dble(nz-1) * (zmax - zmin)
    end do
    f1d = sp%interpolate([nx, ny, nz], [x,y,z])
    f = reshape(f1d, shape(f))

    max_err = 0.d0
    open(15, file="file_3d.dat")
    write(15,*) "x,  y, interpolated, true value"
    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          write(15,"(3f8.4, 2es18.6)") x(i), y(j), z(k), f(i,j,k), f3(x(i),y(j),z(k))
          max_err = max(max_err, abs(f(i,j,k) - f3(x(i),y(j),z(k))))
        end do
      end do
    end do
    close(15)
    write(*,"(a,es18.6)") "Maximum interpolation error: ", max_err
  end subroutine test_3d


  function f1(x)
    real(8), intent(in) :: x
    real(8) :: f1
    f1 = cos(x) * exp(-x**2)
  end function f1

  function f2(x,y)
    real(8), intent(in) :: x, y
    real(8) :: f2
    f2 = cos(x) * cos(y) * exp(-x**2-y**2)
  end function f2

  function f3(x,y,z)
    real(8), intent(in) :: x, y, z
    real(8) :: f3
    f3 = cos(x) * cos(y) * cos(z) * exp(-x**2-y**2-z**2)
  end function f3

  function f4(w,x,y,z)
    real(8), intent(in) :: w, x, y, z
    real(8) :: f4
    f4 = cos(w) * cos(x) * cos(y) * cos(z) * exp(-w**2-x**2-y**2-z**2)
  end function f4
end program test
