module NdSpline
  !
  ! This is for the n-dimensional spline interpolation using the basis functions.
  !
  ! Suppose, we have the function f depending on x, f = f(x), and we know the value of
  ! f(x) at x = x_{i}. Here, i = 1, ..., n.
  ! Then, the interpolated f(x) is given as
  ! f(x) = sum c_{i} phi^{k}_{i}(x), i = 1, ..., n
  ! with the k-th order polynomial basis function phi^{k}_{i}(x) and coefficient c_{i}.
  ! Note that the c_{i} is determined so f(x_{i}) = \sum c_{i} phi^{k}_{i}(x).
  !
  ! With the tensor product form, the generalization to the n dimension is pretty straightforward.
  ! For example, in the 4 dimensional, a function f(w,x,y,z) can be
  ! f(w,x,y,z) = sum c_{ijkl} phi^{kw}_{i}(w) phi^{kx}_{j}(x) phi^{ky}_{k}(y) phi^{kz}_{l}(z)
  !
  ! See G. D. Knott, "Interpolating cubic splines", Springer 2000 for theoretical details.
  !
  ! See the test.F90 in "main" directory how to use.
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: dp = real64
  public ::spline
  private :: abscissa
  private :: abscissa_interpolant
  private :: grid_interpolant
  private :: fin_spline
  private :: init_spline
  private :: interpolate
  private :: interpolate_1d
  private :: set_coefs
  private :: set_coefs_1d
  private :: fin_abscissa
  private :: init_abscissa
  private :: set_knot_vector
  private :: fin_grid_interpolant
  private :: init_grid_interpolant
  private :: fin_abscissa_interpolant
  private :: init_abscissa_interpolant
  private :: find_interval
  private :: basis_function

  type :: abscissa
    real(dp), allocatable :: x(:), t(:)
    integer :: nx = 0
    integer :: nt = 0
    integer :: k  = 0
  contains
    procedure :: init_abscissa
    procedure :: fin_abscissa
    procedure :: set_knot_vector
    procedure :: find_interval
    generic :: init => init_abscissa
    generic :: fin => fin_abscissa
  end type abscissa

  type :: spline
    type(abscissa), allocatable :: Ndrctn(:)
    real(dp), allocatable :: coefs(:)
    integer :: ndim = 0
    integer :: nf_in = 0
  contains
    procedure :: fin_spline
    procedure :: init_spline
    procedure :: set_coefs
    procedure :: interpolate
    generic :: init => init_spline
    generic :: fin => fin_spline
  end type spline

  type :: abscissa_interpolant
    real(dp), allocatable :: x(:)
    integer :: n = 0
  contains
    procedure :: init_abscissa_interpolant
    procedure :: fin_abscissa_interpolant
    generic :: init => init_abscissa_interpolant
    generic :: fin => fin_abscissa_interpolant
  end type abscissa_interpolant

  type :: grid_interpolant
    type(abscissa_interpolant), allocatable :: Ndrctn(:)
    integer :: ndim = 0
    integer :: n_interpolant = 0
  contains
    procedure :: init_grid_interpolant
    procedure :: fin_grid_interpolant
    generic :: init => init_grid_interpolant
    generic :: fin => fin_grid_interpolant
  end type grid_interpolant

contains
  subroutine fin_spline(this)
    class(spline), intent(inout) :: this
    integer :: n
    do n = 1, this%ndim
      call this%Ndrctn(n)%fin()
    end do
    deallocate(this%coefs)
  end subroutine fin_spline

  subroutine init_spline(this, ks, ns, xs, fs)
    ! constructor of spline
    ! inputs:
    ! 'ks' : 1 dimenstionl array for the order of polymial for each abscissa
    ! 'ns' : 1 dimenstionl array for the dimension of each abscissa
    ! 'xs' : 1 dimenstionl array for the points of each abscissa. One can use the concatenated array [x,y,...]
    !           x, y, ... are non dcreasing one dimensional array
    ! 'fs' : 1 dimenstionl array for the points of each abscissa. One can use fs = reshape( f, (nx*ny*...) )
    !         Here f is multidimenstional array f(i,j,...). i,j,... have to be correspond to x,y,..., respectively.
    class(spline), intent(inout) :: this
    integer, intent(in) :: ks(:), ns(:)
    real(dp), intent(in) :: xs(:), fs(:)
    integer :: n, n_start, n_end, n_f

    if(size(ks) /= size(ns)) then
      write(*,"(a)") "Error in constructor of spline: Dimension cannot be detected properly!"
      write(*,"(2i3)") size(ks), size(ns)
      stop
    end if

    n_f = 1
    do n = 1, size(ns)
      n_f = n_f * ns(n)
    end do

    if(size(fs) /= n_f) then
      write(*,"(a)") "Error in constructor of spline: f cannot be detected properly"
      write(*,"(2i8)") n_f, size(fs)
      stop
    end if

    this%ndim = size(ks)
    allocate(this%Ndrctn(this%ndim))
    this%nf_in = n_f
    n_start = 0
    n_end = 0
    do n = 1, this%ndim
      n_start = n_end + 1
      n_end = n_end + ns(n)
      if(ks(n) < 1) then
        write(*,"(a)") "Error in constructor of spline: polynomial order k has to be greater than 0"
        write(*,"(i3)") ks(n)
        stop
      end if
      call this%Ndrctn(n)%init( ks(n), xs(n_start:n_end) )
    end do

    allocate(this%coefs(size(fs)))
    this%coefs = fs
    call this%set_coefs()
  end subroutine init_spline

  function interpolate(this, ns, xs) result(f)
    ! interpolation using the matrix-matrix multiplication
    ! inputs:
    ! 'ns': 1 dimenstional array specifying the dimenstion
    ! 'xs': 1 dimenstional array specifying the interpolants
    ! output:
    ! 'f': 1 dimenstional array
    class(spline), intent(inout) :: this
    integer, intent(in) :: ns(:)
    real(dp), intent(in) :: xs(:)
    real(dp), allocatable :: f(:)
    type(grid_interpolant) :: grid
    real(dp), allocatable :: tmp_org(:,:), tmp_int(:,:)
    integer :: i, n

    call grid%init(ns, xs)
    if(this%ndim /= grid%ndim) then
      write(*,"(a)") "Error in interpolate of spline: Dimension would be wrong!"
      return
    end if

    n = this%nf_in
    do i = 1, this%ndim
      allocate(tmp_org(this%Ndrctn(i)%nx, n / this%Ndrctn(i)%nx))

      if(i == 1) then
        tmp_org = reshape(this%coefs, shape(tmp_org))
      end if

      if(i /= 1) then
        tmp_org = reshape(tmp_int, shape(tmp_org))
        deallocate(tmp_int)
      end if
      allocate(tmp_int(n / this%Ndrctn(i)%nx, grid%Ndrctn(i)%n))

      call interpolate_1d(this%Ndrctn(i), grid%Ndrctn(i)%x, tmp_org, tmp_int)
      deallocate(tmp_org)
      n = n * grid%Ndrctn(i)%n / this%Ndrctn(i)%nx
    end do
    allocate(f(grid%n_interpolant))
    f = reshape(tmp_int, shape(f))
    deallocate(tmp_int)
    call grid%fin()
  end function interpolate

  subroutine interpolate_1d(this, x, coefs, f)
    !
    !  f_{n} = \sum c_{i} phi_{ni}
    !
    type(abscissa), intent(in) :: this
    real(dp), intent(in) :: x(:), coefs(:,:)
    real(dp), intent(inout) :: f(:,:)
    real(dp), allocatable :: Bmat(:,:)
    integer :: n, m, i, j
    n = size(this%x)
    m = size(x)
    allocate(BMat(m,n))
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, m
        Bmat(j,i) = basis_function(this, i, this%k, x(j))
      end do
    end do
    !$omp end do
    !$omp end parallel
    call dgemm("t", "t", size(coefs,2), m, n, 1.0_dp, coefs, n, BMat, m, 0.0_dp, f, size(coefs,2))
  end subroutine interpolate_1d

  subroutine set_coefs(this)
    !
    ! Here, the coefficient c is determine.
    !
    class(spline), intent(inout) :: this
    real(dp), allocatable :: tmp_coef(:,:), tmp(:,:)
    integer :: i

    do i = 1, this%ndim

      allocate(tmp_coef(this%Ndrctn(i)%nx, this%nf_in / this%Ndrctn(i)%nx))
      if(i == 1) then
        tmp_coef = reshape(this%coefs, shape(tmp_coef))
      end if

      if(i /= 1) then
        tmp_coef = reshape(transpose(tmp), shape(tmp_coef))
        deallocate(tmp)
      end if
      allocate(tmp(this%Ndrctn(i)%nx, this%nf_in / this%Ndrctn(i)%nx))

      call set_coefs_1d(this%Ndrctn(i), tmp_coef, tmp)
      deallocate(tmp_coef)
    end do
    this%coefs = reshape(transpose(tmp), shape(this%coefs))
  end subroutine set_coefs

  subroutine set_coefs_1d(this, f, coefs)
    !
    ! determine c_{i} so f_{n} = \sum c_{i} phi_{ni}.
    ! Here f_{n} indicates the initially given point.
    !
    type(abscissa), intent(in) :: this
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(inout) :: coefs(:,:)
    integer :: i, j, info
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: Bmat(:,:)

    allocate(BMat(this%nx, this%nx))
    allocate(ipiv(this%nx))
    !$omp parallel
    !$omp do private(i,j)
    do i = 1, this%nx
      do j = 1, this%nx
        BMat(j,i) = basis_function(this, i, this%k, this%x(j))
      end do
    end do
    !$omp end do
    !$omp end parallel
    coefs = f
    call dgesv(this%nx, size(coefs,2), BMat, this%nx, ipiv, coefs, this%nx, info)
    deallocate(BMat, ipiv)
  end subroutine set_coefs_1d

  subroutine fin_abscissa(this)
    class(abscissa), intent(inout) :: this
    deallocate(this%x)
    deallocate(this%t)
  end subroutine fin_abscissa

  subroutine init_abscissa(this, k, x)
    class(abscissa), intent(inout) :: this
    integer, intent(in) :: k
    real(dp), intent(in) :: x(:)

    this%k = k
    this%nx = size(x)
    this%nt = this%k + this%nx
    allocate(this%x(this%nx))
    allocate(this%t(this%nt))
    this%x = x
    call this%set_knot_vector()
  end subroutine init_abscissa

  subroutine set_knot_vector(this)
    class(abscissa), intent(inout) :: this
    integer :: n, k, i
    real(dp), allocatable :: tmp(:)

    n = this%nx
    k = this%k
    this%t(:) = 0.d0
    this%t(:k) = this%x(1)
    this%t(n+1:) = this%x(n)

    if(mod(k,2) == 0) then
      this%t(k+1:n) = this%x(k/2+1:n-k/2)
      return
    end if

    if(mod(k,2) == 1) then
      allocate(tmp( n-k ))
      do i = 1, n-k
        tmp(i) = 0.5d0 * ( this%x(k/2+i) + this%x(k/2+i+1) )
      end do
      this%t(k+1:n) = tmp(:)
      return
    end if
  end subroutine set_knot_vector

  subroutine fin_grid_interpolant(this)
    class(grid_interpolant), intent(inout) :: this
    integer :: i
    do i = 1, this%ndim
      call this%Ndrctn(i)%fin()
    end do
    this%n_interpolant = 0
    this%ndim = 0
  end subroutine fin_grid_interpolant

  subroutine init_grid_interpolant(this, ns, xs)
    class(grid_interpolant), intent(inout) :: this
    integer, intent(in) :: ns(:)
    real(dp), intent(in) :: xs(:)
    integer :: i, n_start, n_end

    this%ndim = size(ns)
    allocate(this%Ndrctn(this%ndim))
    n_start = 0
    n_end = 0
    do i = 1, this%ndim
      n_start = n_end + 1
      n_end = n_end + ns(i)
      call this%Ndrctn(i)%init( xs(n_start:n_end) )
    end do
    this%n_interpolant = 1
    do i = 1, this%ndim
      this%n_interpolant = this%n_interpolant * this%Ndrctn(i)%n
    end do
  end subroutine init_grid_interpolant

  subroutine fin_abscissa_interpolant(this)
    class(abscissa_interpolant), intent(inout) :: this
    deallocate(this%x)
    this%n = 0
  end subroutine fin_abscissa_interpolant

  subroutine init_abscissa_interpolant(this, x)
    class(abscissa_interpolant), intent(inout) :: this
    real(dp), intent(in) :: x(:)

    this%n = size(x)
    allocate(this%x(this%n))
    this%x = x
  end subroutine init_abscissa_interpolant

  function find_interval(this, x) result(idx)
    ! return the index i so t(i) <= x <= t(i+1)
    ! Note that x has to satisfy t(1) <= x <= t(n+k)
    class(abscissa), intent(in) :: this
    real(dp), intent(in) :: x
    integer :: idx, i, k, nt

    k = this%k
    nt = size(this%t)
    idx = 0
    if(x < minval(this%t) .or. x > maxval(this%t)) return

    if(this%t(k) <= x .and. x <= this%t(k+1)) then
      idx = k
      return
    end if

    if(this%t(nt-k) <= x .and. x <= this%t(nt-k+1)) then
      idx = nt-k
      return
    end if

    do i = 1, nt - 1
      if( this%t(i) <= x .and. x <= this%t(i+1) ) then
        idx = i
        return
      end if
    end do
  end function find_interval

  recursive function basis_function(this, i, k, x) result(f)
    !
    ! basis function phi^{k}_{i}(x) using the Cox-de Boor recursion formula
    !
    type(abscissa), intent(in) :: this
    integer, intent(in) :: i, k
    real(dp), intent(in) :: x
    real(dp) :: f, c1, c2

    f = 0.0_dp
    if( i + k > size(this%t)) return
    if( k == 1 ) then
      if(i == this%find_interval(x)) f = 1.0_dp
      return
    end if
    c1 = (x - this%t(i)) / (this%t(i+k-1) - this%t(i))
    c2 = (this%t(i+k)-x) / (this%t(i+k) - this%t(i+1))
    if( abs(this%t(i+k-1) - this%t(i)) < 1.d-8 ) c1 = 0.0_dp
    if( abs(this%t(i+k) - this%t(i+1)) < 1.d-8 ) c2 = 0.0_dp
    if( abs(this%t(i+k-1) - this%t(i)) < 1.d-8 .and. abs(x - this%t(i)) < 1.d-8 ) c1 = 1.0_dp
    if( abs(this%t(i+k) - this%t(i+1)) < 1.d-8 .and. abs(x - this%t(i+k)) < 1.d-8 ) c2 = 1.0_dp
    f = basis_function(this, i, k-1, x) * c1 + basis_function(this, i+1, k-1, x) * c2
  end function basis_function
end module NdSpline


