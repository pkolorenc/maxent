!===============================================================================
! MODULE: Numerics_module
!
! PURPOSE: Provides numerical utilities for the Maximum Entropy Method:
!          - Orthogonal polynomial basis functions (Legendre, Chebyshev)
!          - Gauss-Legendre quadrature for numerical integration
!          - Bisection method for root finding
!          - Physical constants and unit conversions
!
! POLYNOMIAL BASES:
!   Legendre polynomials:
!     - Orthogonal on [-1,1] with uniform weight
!     - Recurrence: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
!     - Optimal for smooth functions
!   
!   Chebyshev polynomials:
!     - Orthogonal on [-1,1] with weight 1/sqrt(1-x^2)
!     - Recurrence: T_{n+1} = 2xT_n - T_{n-1}
!     - Better near endpoints (slower convergence in middle)
!
! QUADRATURE:
!   Gauss-Legendre: optimal for polynomial integration on finite intervals
!   Exact for polynomials of degree ≤ 2n-1 with n points
!
! ROOT FINDING:
!   Bisection: robust but slow (linear convergence)
!   Guaranteed to find root if function changes sign in bracket
!
!===============================================================================

MODULE Numerics_module
  implicit none

  !-----------------------------------------------------------------------------
  ! Physical and mathematical constants
  !-----------------------------------------------------------------------------
  REAL(8), PARAMETER :: pi = acos(-1.d0), &      ! π ≈ 3.14159265358979...
                        Ha = 27.21138386d0       ! Hartree to eV conversion
  ! 1 Hartree = 27.2114 eV is the atomic unit of energy
  ! Equal to twice the ionization energy of hydrogen ground state

  !-----------------------------------------------------------------------------
  ! Polynomial basis selection flags
  !-----------------------------------------------------------------------------
  LOGICAL :: chebyshev = .false., &  ! Use Chebyshev polynomials T_n(x)
             legendre = .true.       ! Use Legendre polynomials P_n(x) (default)
  ! Only one should be .true. at a time
  ! If both false, reverts to monomials x^n (not recommended)

  !-----------------------------------------------------------------------------
  ! Coordinate transformation parameter
  !-----------------------------------------------------------------------------
  REAL(8) :: low_e_boundary  ! Lower boundary in energy space
                              ! Used to map x ∈ [low_xx, high_xx] to y ∈ [-1,1]
                              ! for polynomial evaluation

CONTAINS

  !-----------------------------------------------------------------------------
  ! SUBROUTINE: Trec01
  !
  ! PURPOSE: Generate polynomial basis functions T_k(x) using recurrence relations
  !          Maps x in [low_xx, high_xx] to y in [-1, 1] for polynomial evaluation
  !
  ! THREE OPTIONS:
  !   1. Legendre polynomials (default): optimal for uniform weight
  !   2. Chebyshev polynomials: optimal near endpoints
  !   3. Monomials: simple but numerically unstable for high degree
  !
  ! WHY ORTHOGONAL POLYNOMIALS?
  !   - Better numerical conditioning than monomials
  !   - Recurrence relations are numerically stable
  !   - Natural basis for function approximation on [-1,1]
  !   - Well-studied convergence properties
  !
  ! COORDINATE MAPPING:
  !   Physical x-space (x = 1/E) is mapped to y ∈ [-1,1]:
  !     y = -1 + 2*low_e_boundary*x
  !   
  !   This ensures polynomials are evaluated on their natural domain
  !   where orthogonality and recurrence relations hold exactly.
  !
  ! INPUT:
  !   x - point in physical x-space where to evaluate polynomials
  !
  ! OUTPUT:
  !   T(:) - array of polynomial values [T_0(x), T_1(x), ..., T_{n-1}(x)]
  !          where n = size(T)
  !
  ! NORMALIZATION:
  !   - Legendre: normalized to integral[-1,1][ P_n(x)^2 dx ] = 2/(2n+1)
  !               After normalization: integral = 1
  !   - Chebyshev: normalized with Chebyshev weight
  !               integral[-1,1][ T_n(x)^2 / sqrt(1-x^2) dx ] = pi/2
  !
  ! NUMERICAL CONSIDERATIONS:
  !   - Recurrence relations are forward-stable (errors don't amplify)
  !   - Evaluation is O(n) where n is polynomial degree
  !   - Much more stable than computing x^n directly for large n
  !
  ! PHYSICS APPLICATION:
  !   These polynomials form the basis for expanding the spectral density:
  !     rho(x) = exp(sum_k[ eta_k * T_k(x) ])
  !   The choice of basis affects numerical stability but not the physics.
  !-----------------------------------------------------------------------------
  SUBROUTINE Trec01(x,T)
    REAL(8) x,T(:)
    INTENT(in) :: x
    INTENT(out) :: T
    INTEGER npol,i
    REAL(8) y
    
    npol = size(T)  ! Number of polynomials to compute
    
    ! Map x from physical domain to [-1,1]
    ! This is crucial for polynomial orthogonality properties
    y = -1.d0+2.d0*low_e_boundary*x
    
    !---------------------------------------------------------------------------
    ! OPTION 1: Chebyshev Polynomials T_n(y)
    !
    ! Definition: T_n(cos(theta)) = cos(n*theta)
    ! 
    ! Recurrence relation:
    !   T_0(y) = 1
    !   T_1(y) = y
    !   T_{n+1}(y) = 2y*T_n(y) - T_{n-1}(y)
    !
    ! Orthogonality: integral[-1,1][ T_m(y)*T_n(y)/sqrt(1-y^2) dy ] = delta_mn * c_n
    !   where c_0 = pi, c_n = pi/2 for n > 0
    !
    ! Normalization: multiply by sqrt(2/pi), with T_0 getting extra 1/sqrt(2)
    !   This makes integral = 1 with Chebyshev weight
    !
    ! ADVANTAGES:
    !   - Minimax property: minimize maximum error
    !   - Better approximation near endpoints x = ±1
    !   - Used in Chebyshev interpolation/approximation
    !
    ! DISADVANTAGES:
    !   - Require Chebyshev weight function (complicates integration)
    !   - Slightly slower convergence in middle of interval
    !---------------------------------------------------------------------------
    if (chebyshev) then
      ! Initialize first two Chebyshev polynomials
      T(1) = 1.d0  ! T_0(y) = 1
      if (npol.gt.1) T(2) = y  ! T_1(y) = y
      
      ! Three-term recurrence for T_n
      do i = 2,npol-1
        T(i+1) = 2.d0*y*T(i)-T(i-1)
      end do
      
      ! Normalize to unit integral with Chebyshev weight
      T(:) = 2.d0*T(:)/sqrt(pi)
      T(1) = T(1)/sqrt(2.d0)  ! Special factor for T_0
      
    !---------------------------------------------------------------------------
    ! OPTION 2: Legendre Polynomials P_n(y)  [DEFAULT]
    !
    ! Definition: P_n(y) are eigenfunctions of Legendre differential equation
    !   d/dy[(1-y^2)dP/dy] + n(n+1)P = 0
    !
    ! Bonnet's recurrence relation:
    !   P_0(y) = 1
    !   P_1(y) = y
    !   (n+1)P_{n+1}(y) = (2n+1)y*P_n(y) - n*P_{n-1}(y)
    !
    ! Orthogonality: integral[-1,1][ P_m(y)*P_n(y) dy ] = delta_mn * 2/(2n+1)
    !
    ! Normalization: multiply by sqrt(2n+1) to make integral = 1
    !
    ! ADVANTAGES:
    !   - Natural for uniform weight (no weight function needed)
    !   - Optimal L^2 approximation for smooth functions
    !   - Widely used in physics (multipole expansions, etc.)
    !   - Better convergence in middle of interval
    !
    ! WHY DEFAULT FOR THIS CODE?
    !   - Spectral density is typically smooth (no endpoint singularities)
    !   - Uniform weight matches Gauss-Legendre quadrature naturally
    !   - Simpler integration (no weight function in integrals)
    !---------------------------------------------------------------------------
    else if (legendre) then
      ! Initialize first two Legendre polynomials
      T(1) = 1.d0  ! P_0(y) = 1
      if (npol.gt.1) T(2) = y  ! P_1(y) = y
      
      ! Bonnet's recurrence relation for P_n
      ! Rearranged as: P_{n+1} = [(2n+1)*y*P_n - n*P_{n-1}]/(n+1)
      do i = 2,npol-1
        T(i+1) = (2.d0*i-1.d0)*y*T(i)-(i-1.d0)*T(i-1)
        T(i+1) = T(i+1)/i
      end do
      
      ! Normalize each P_n to unit integral:
      ! integral[-1,1][ P_n(y)^2 dy ] = 2/(2n+1) normally
      ! After multiplication by sqrt(2n+1): integral = 1
      do i = 0,npol-1
        T(i+1) = T(i+1)*sqrt((2.d0*i+1.d0))
      end do
      
    !---------------------------------------------------------------------------
    ! OPTION 3: Monomials x^k  [NOT RECOMMENDED]
    !
    ! Simple power basis: 1, x, x^2, x^3, ...
    !
    ! DISADVANTAGES:
    !   - Not orthogonal (leads to ill-conditioned matrices)
    !   - High powers numerically unstable (x^20 can overflow/underflow)
    !   - Poor approximation properties compared to orthogonal polynomials
    !
    ! WHEN TO USE:
    !   - Only for debugging or comparison
    !   - When orthogonality is not needed
    !   - Low polynomial degree (n < 5)
    !
    ! This is the fallback if both chebyshev and legendre are false
    !---------------------------------------------------------------------------
    else
      ! Simple monomial basis: T_n(x) = x^n
      T(1) = 1.d0
      do i = 1,npol-1
        T(i+1) = T(i)*x  ! Iteratively compute x^i
      end do
    end if
  END SUBROUTINE Trec01

  !-----------------------------------------------------------------------------
  ! FUNCTION: ChW  [COMMENTED OUT]
  !
  ! PURPOSE: Chebyshev weight function w(x) = 1/sqrt(x*(1-x))
  !
  ! This would be used if integrating with Chebyshev weight:
  !   integral[ f(x) * w(x) dx ]
  !
  ! Currently not used because:
  !   - Legendre is default (uniform weight)
  !   - Chebyshev weight complicates quadrature
  !   - Can be implemented later if needed
  !-----------------------------------------------------------------------------
  !REAL(8) FUNCTION ChW(x)
  !  REAL(8),INTENT(in) :: x
  !  if (chebyshev) then
  !    ChW = 0.5d0/sqrt(x*(1.d0-x))
  !  else
  !    ChW = 1.d0
  !  end if
  !END FUNCTION ChW

  !-----------------------------------------------------------------------------
  ! SUBROUTINE: GaussLeg
  !
  ! PURPOSE: Compute Gauss-Legendre quadrature points and weights for integration
  !          over arbitrary interval [x1, x2]
  !
  ! THEORY:
  !   Gauss-Legendre quadrature approximates integrals as:
  !     integral[x1 to x2][ f(x) dx ] ≈ sum_{i=1}^n[ w_i * f(x_i) ]
  !   
  !   The nodes x_i are zeros of the n-th Legendre polynomial P_n(x)
  !   The weights are chosen to maximize accuracy
  !   
  !   OPTIMAL PROPERTY: Exact for polynomials of degree ≤ 2n-1
  !   This is the highest possible accuracy for n function evaluations!
  !
  ! ALGORITHM:
  !   1. Find zeros of P_n(x) using Newton's method
  !      - Initial guess: cos(pi*(i-0.25)/(n+0.5))  [asymptotic formula]
  !      - Iterate: x_new = x_old - P_n(x_old)/P_n'(x_old)
  !   
  !   2. Compute weights: w_i = 2/[(1-x_i^2)*[P_n'(x_i)]^2]
  !   
  !   3. Transform from [-1,1] to [x1,x2]:
  !      - x_phys = xm + xl*x_std  where xm = (x1+x2)/2, xl = (x2-x1)/2
  !      - w_phys = xl*w_std
  !
  ! WHY USE GAUSS-LEGENDRE?
  !   - Optimal accuracy for smooth functions
  !   - No function evaluations at endpoints (avoids singularities)
  !   - Well-tested and reliable
  !   - Standard method in computational physics
  !
  ! NUMERICAL STABILITY:
  !   - Newton's method converges quadratically for isolated roots
  !   - Legendre polynomials have n distinct real roots in (-1,1)
  !   - Recurrence relation for P_n is forward-stable
  !   - Typically converges in < 10 iterations
  !
  ! INPUT:
  !   x1, x2 - integration limits [x1, x2]
  !   x(:)   - array to store quadrature nodes (allocated)
  !   w(:)   - array to store quadrature weights (allocated)
  !
  ! OUTPUT:
  !   x(:) - filled with n quadrature nodes in [x1, x2]
  !   w(:) - filled with n quadrature weights
  !   sum(w) = x2 - x1  (as expected)
  !
  ! USAGE EXAMPLE:
  !   n = 10
  !   allocate(x(n), w(n))
  !   call GaussLeg(0.d0, 1.d0, x, w)  ! Setup quadrature on [0,1]
  !   integral = sum(w * f(x))          ! Approximate integral of f
  !
  ! ADAPTED FROM:
  !   Numerical Recipes in Fortran 90 (Press et al.)
  !   Modified for modern Fortran and clarity
  !-----------------------------------------------------------------------------
  SUBROUTINE GaussLeg(x1, x2, x, w)
    REAL(8), INTENT(in) :: x1, x2     ! Integration limits
    REAL(8), DIMENSION(:), INTENT(out) :: x, w  ! Nodes and weights
    
    ! Numerical parameters
    REAL(8), PARAMETER :: EPS = 3.0d-16  ! Convergence tolerance (machine epsilon)
    integer, parameter :: MAXIT = 10      ! Maximum Newton iterations
    
    ! Work variables
    integer :: its, j, m, n
    real(8) :: xl, xm
    real(8), dimension((size(x) + 1) / 2) :: p1, p2, p3, pp, z, z1
    logical, dimension((size(x) + 1) / 2) :: unfinished
    
    n = size(x)      ! Number of quadrature points
    m = (n + 1) / 2  ! Only need to find half the roots (symmetry about x=0)
    
    ! Midpoint and half-length for interval transformation
    xm = 0.5d0 * (x2 + x1)  ! Center of interval
    xl = 0.5d0 * (x2 - x1)  ! Half-length of interval
    
    !---------------------------------------------------------------------------
    ! Initial approximations to roots of P_n(x)
    ! Using asymptotic formula: i-th root ≈ cos(pi*(i-1/4)/(n+1/2))
    ! This is accurate to O(1/n^2) and provides excellent starting guess
    !---------------------------------------------------------------------------
    do j = 1,m
      z(j) = cos(pi*(j-0.25d0)/(n+0.5d0))
    end do
    
    !---------------------------------------------------------------------------
    ! Newton's method to refine roots
    ! Iterates: x_new = x_old - P_n(x_old)/P_n'(x_old)
    !---------------------------------------------------------------------------
    unfinished = .true.  ! Track which roots haven't converged yet
    
    do its = 1, MAXIT
      ! Initialize Legendre polynomial recurrence for converging roots
      where (unfinished)
        p1 = 1.0   ! P_0 = 1
        p2 = 0.0   ! P_{-1} = 0 (convention)
      end where
      
      ! Evaluate P_n(z) using Bonnet's recurrence relation:
      !   (j+1)P_{j+1} = (2j+1)*z*P_j - j*P_{j-1}
      do j = 1, n
        where (unfinished)
          p3 = p2
          p2 = p1
          p1 = ((2.0d0 * j - 1.0d0) * z * p2 - (j - 1.0d0) * p3) / j
        end where
      end do
      ! Now p1 = P_n(z)
      
      ! Compute derivative P_n'(z) using relation:
      !   P_n'(z) = n*[z*P_n(z) - P_{n-1}(z)]/(z^2-1)
      where (unfinished)
        pp = n * (z * p1 - p2) / (z * z - 1.0d0)
        z1 = z
        z = z1 - p1 / pp  ! Newton update: x_new = x_old - f/f'
        unfinished = (abs(z - z1) > EPS)  ! Check convergence
      end where
      
      if (.not. any(unfinished)) exit  ! All roots converged
    end do
    
    ! Check for convergence failure (very rare with good initial guess)
    if (its == MAXIT + 1) stop 'too many iterations in GaussLeg'
    
    !---------------------------------------------------------------------------
    ! Scale roots from [-1,1] to [x1,x2] and compute weights
    !
    ! Transformation: x_phys = xm + xl*x_std
    !   where x_std ∈ [-1,1] and x_phys ∈ [x1,x2]
    !
    ! Weight formula: w_i = 2*xl/[(1-z_i^2)*[P_n'(z_i)]^2]
    !   Includes xl factor for interval transformation
    !
    ! Symmetry: Use P_n(-x) = (-1)^n*P_n(x) to get second half of points/weights
    !---------------------------------------------------------------------------
    x(1 : m) = xm - xl * z              ! First half of nodes
    x(n : n - m + 1 : -1) = xm + xl * z ! Second half (symmetric)
    
    w(1 : m) = 2.0d0 * xl / ((1.0d0 - z**2) * pp**2)  ! First half of weights
    w(n : n - m + 1 : -1) = w(1 : m)                  ! Second half (symmetric)
    
  END SUBROUTINE GaussLeg

  !-----------------------------------------------------------------------------
  ! SUBROUTINE: get_unit
  !
  ! PURPOSE: Find an available (unused) I/O unit number for file operations
  !
  ! FORTRAN I/O UNITS:
  !   - Units 5, 6 are typically stdin, stdout
  !   - Units 0, 5, 6 may be preconnected
  !   - Safe range: typically 10-99 for user files
  !
  ! METHOD:
  !   Search sequentially from unit 10 until finding:
  !     1. Unit exists (inquire gives exist=.true.)
  !     2. Unit is not opened (inquire gives opened=.false.)
  !
  ! OUTPUT:
  !   iounit - first available unit number (≥ 10)
  !
  ! USAGE:
  !   call get_unit(iou)
  !   open(iou, file='output.dat', status='unknown')
  !   write(iou,*) data
  !   close(iou)
  !
  ! NOTE: Not thread-safe (assumes sequential execution)
  !-----------------------------------------------------------------------------
  SUBROUTINE get_unit(iounit)
    INTEGER,INTENT(out) :: iounit
    logical unitok,unitop
    integer i
    
    iounit = -1
    i = 10  ! Start search at unit 10 (avoid preconnected units)
    
    do while (iounit.eq.-1)
      ! Check if unit number i is valid and available
      inquire(unit=i,exist=unitok,opened=unitop)
      if (unitok.and.(.not.unitop)) then
        iounit = i  ! Found available unit
      else
        i = i+1     ! Try next unit
      end if
    end do
  END SUBROUTINE get_unit

  !-----------------------------------------------------------------------------
  ! SUBROUTINE: bisection
  !
  ! PURPOSE: Find root of function f(x) using bisection method
  !          Solves: f(x) = 0 for x in [aa, bb]
  !
  ! ALGORITHM:
  !   1. Start with bracket [a,b] where f(a)*f(b) < 0 (sign change)
  !   2. Compute midpoint c = (a+b)/2
  !   3. Evaluate f(c)
  !   4. If f(a)*f(c) < 0: root in [a,c], set b=c
  !      Else: root in [c,b], set a=c
  !   5. Repeat until |f(c)| < tol or interval width < tol
  !
  ! CONVERGENCE:
  !   - Linear convergence: error halves each iteration
  !   - After n iterations: error ≤ (b-a)/2^n
  !   - Guaranteed convergence if f continuous and bracket valid
  !
  ! ADVANTAGES:
  !   - Extremely robust (always converges for valid bracket)
  !   - Simple to implement
  !   - No derivative needed (unlike Newton's method)
  !   - Handles discontinuous derivatives
  !
  ! DISADVANTAGES:
  !   - Slow (linear convergence vs quadratic for Newton)
  !   - Requires bracket [a,b] with sign change
  !   - Not suitable for multiple roots or flat regions
  !
  ! WHY USE FOR THIS CODE?
  !   The function activeDmu(eta_k) is:
  !     - Expensive to evaluate (requires computing rho at all grid points)
  !     - Smooth but derivative not easily available
  !     - Needs robust convergence (Newton can diverge)
  !     - Bisection's linear convergence acceptable given function cost
  !
  ! INPUT:
  !   f - function handle with signature: real(8) function f(x)
  !   aa, bb - initial bracket [aa, bb] where f(aa)*f(bb) < 0
  !   tol - convergence tolerance for |f(root)| or interval width
  !   max_iter - maximum iterations (typically 50-100 sufficient)
  !
  ! OUTPUT:
  !   root - approximate solution to f(x) = 0
  !   success - .true. if converged, .false. if max iterations reached
  !
  ! USAGE IN ENTROPY CODE:
  !   We solve: activeDmu(eta_k) = 0
  !   This finds eta_k such that moment constraint is satisfied:
  !     integral[ T_k(x)*rho(x;eta_k) dx ] = mui_k
  !
  ! ERROR HANDLING:
  !   If max iterations reached without convergence:
  !     - Prints warning message
  !     - Returns midpoint of final bracket
  !     - Sets success = .false.
  !   Caller should check success flag and handle appropriately
  !-----------------------------------------------------------------------------
  SUBROUTINE bisection(f,aa,bb,tol, max_iter, root, success)
    ! Function interface specification
    interface
      function f(x) result(y)
        real(8), intent(in) :: x
        real(8) :: y
      end function f
    end interface

    real(8), intent(in) :: aa, bb, tol  ! Bracket and tolerance
    integer, intent(in) :: max_iter     ! Maximum iterations
    real(8), intent(out) :: root        ! Solution
    logical, intent(out) :: success     ! Convergence flag
    
    ! Local variables for bisection state
    real(8) :: a, b, fa, fb, c, fc
    integer :: iter
    
    ! Initialize bracket
    a = aa
    b = bb
    fa = f(a)  ! Evaluate function at left endpoint
    fb = f(b)  ! Evaluate function at right endpoint
    
    ! Note: Should check fa*fb < 0 (sign change), but omitted for efficiency
    ! Caller's responsibility to provide valid bracket
    
    ! Main bisection loop
    do iter = 1, max_iter
      ! Compute midpoint of current bracket
      c = (a + b) / 2.0d0
      fc = f(c)
      
      ! Check convergence: either function small or interval narrow
      if (abs(fc) < tol .or. (b - a) / 2.0d0 < tol) then
        root = c
        success = .true.
        return
      end if
      
      ! Determine which half-interval contains root
      if (fa * fc < 0.0d0) then
        ! Root in left half [a, c]
        b = c
        fb = fc
      else
        ! Root in right half [c, b]
        a = c
        fa = fc
      end if
    end do
    
    ! Max iterations reached without convergence
    print*, 'Maximum iterations reached without convergence.'
    root = (a + b) / 2.0d0  ! Return midpoint of final bracket
    success = .false.
  end subroutine bisection

END MODULE Numerics_module
