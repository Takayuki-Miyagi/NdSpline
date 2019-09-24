program test
  use omp_lib
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
    call test_4d()
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
    real(8) :: ti

    n = 30
    allocate(x(n), f(n))
    do i = 1, n
      x(i) = xmin + dble(i-1) / dble(n-1) * (xmax - xmin)
    end do

    do i = 1, n
      f(i) = f1(x(i))
    end do

    ti = omp_get_wtime()
    call sp%init([k],[n],x,f)
    write(*,"(a,f12.6,a)") "constructor: ", omp_get_wtime() - ti, " sec"
    deallocate(x, f)

    n = 100
    allocate(x(n))
    do i = 1, n
      x(i) = xmin + dble(i-1) / dble(n-1) * (xmax - xmin)
    end do

    ti = omp_get_wtime()
    f = sp%interpolate([n], x)
    write(*,"(a,f12.6,a)") "interpolate: ", omp_get_wtime() - ti, " sec"

    call sp%fin()

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
    real(8) :: ti

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
    ti = omp_get_wtime()
    call sp%init([kx,ky],[nx,ny],[x,y],f1d)
    write(*,"(a,f12.6,a)") "constructor: ", omp_get_wtime() - ti, " sec"
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
    ti = omp_get_wtime()
    f = reshape(f1d, shape(f))
    write(*,"(a,f12.6,a)") "interpolate: ", omp_get_wtime() - ti, " sec"

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
    real(8) :: xmin = 0.d0, xmax = 8.d0
    real(8) :: ymin = 0.d0, ymax = 8.d0
    real(8) :: zmin = 0.d0, zmax = 8.d0
    integer :: nx, ny, nz, kx=4, ky=4, kz = 4
    real(8), allocatable :: x(:), y(:), z(:), f(:,:,:), f1d(:)
    real(8) :: max_err
    integer :: i, j, k
    real(8) :: ti

    nx = 15
    ny = 15
    nz = 15
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
    ti = omp_get_wtime()
    call sp%init([kx,ky,kz],[nx,ny,nz],[x,y,z],f1d)
    write(*,"(a,f12.6,a)") "constructor: ", omp_get_wtime() - ti, " sec"
    deallocate(x, y, z, f, f1d)

    nx = 60
    ny = 60
    nz = 60
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
    ti = omp_get_wtime()
    f1d = sp%interpolate([nx, ny, nz], [x,y,z])
    write(*,"(a,f12.6,a)") "interpolate: ", omp_get_wtime() - ti, " sec"
    f = reshape(f1d, shape(f))

    max_err = 0.d0
    open(15, file="file_3d.dat")
    write(15,*) "x,  y, z,  interpolated, true value"
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

  subroutine test_4d()
    type(spline) :: sp
    real(8) :: wmin = 0.d0, wmax = 8.d0
    real(8) :: xmin = 0.d0, xmax = 8.d0
    real(8) :: ymin = 0.d0, ymax = 8.d0
    real(8) :: zmin = 0.d0, zmax = 8.d0
    integer :: nw, nx, ny, nz, kw=4, kx=4, ky=4, kz = 4
    real(8), allocatable :: w(:), x(:), y(:), z(:), f(:,:,:,:), f1d(:)
    real(8) :: max_err
    integer :: i, j, k, l
    real(8) :: ti

    nw = 15
    nx = 15
    ny = 15
    nz = 15
    allocate(w(nw), x(nx), y(ny), z(nz), f(nw,nx,ny,nz), f1d(nw*nx*ny*nz))
    do i = 1, nw
      w(i) = wmin + dble(i-1) / dble(nw-1) * (wmax - wmin)
    end do
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do
    do i = 1, nz
      z(i) = zmin + dble(i-1) / dble(nz-1) * (zmax - zmin)
    end do

    do i = 1, nw
      do j = 1, nx
        do k = 1, ny
          do l = 1, nz
            f(i,j,k,l) = f4(w(i),x(j),y(k),z(l))
          end do
        end do
      end do
    end do

    f1d = reshape(f, shape(f1d))
    ti = omp_get_wtime()
    call sp%init([kw,kx,ky,kz],[nw,nx,ny,nz],[w,x,y,z],f1d)
    write(*,"(a,f12.6,a)") "constructor: ", omp_get_wtime() - ti, " sec"
    deallocate(w, x, y, z, f, f1d)

    nw = 60
    nx = 60
    ny = 60
    nz = 60
    allocate(w(nw),x(nx),y(ny),z(nz),f(nw,nx,ny,nz))
    do i = 1, nw
      w(i) = wmin + dble(i-1) / dble(nw-1) * (wmax - wmin)
    end do
    do i = 1, nx
      x(i) = xmin + dble(i-1) / dble(nx-1) * (xmax - xmin)
    end do
    do i = 1, ny
      y(i) = ymin + dble(i-1) / dble(ny-1) * (ymax - ymin)
    end do
    do i = 1, nz
      z(i) = zmin + dble(i-1) / dble(nz-1) * (zmax - zmin)
    end do
    ti = omp_get_wtime()
    f1d = sp%interpolate([nw, nx, ny, nz], [w, x, y, z])
    write(*,"(a,f12.6,a)") "interpolate: ", omp_get_wtime() - ti, " sec"
    f = reshape(f1d, shape(f))

    max_err = 0.d0
    open(15, file="file_4d.dat")
    write(15,*) "w,  x,  y, z,  interpolated, true value"
    do i = 1, nw
      do j = 1, nx
        do k = 1, ny
          do l = 1, nz
            !write(15,"(4f8.4, 2es18.6)") w(i), x(j), y(k), z(l), f(i,j,k,l), f4(w(i),x(j),y(k),z(l))
            max_err = max(max_err, abs(f(i,j,k,l) - f4(w(i),x(j),y(k),z(l))))
          end do
        end do
      end do
    end do
    close(15)
    write(*,"(a,es18.6)") "Maximum interpolation error: ", max_err
  end subroutine test_4d


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
