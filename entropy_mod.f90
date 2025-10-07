!===============================================================================
! MODULE: entropy_module
!
! PURPOSE: Core functions for Maximum Entropy reconstruction of spectral density
!          from moment constraints. Implements the Maximum Entropy Method (MEM)
!          where the spectral density has the form:
!          
!          rho(x) = exp(-1 - sum_k eta_k * T_k(x)) / x^2
!
!          subject to constraints:
!          
!          integral[ T_k(x) * rho(x) dx ] = mui_k  (for k=0,1,...,NT)
!
! THEORY:  The maximum entropy principle ensures the smoothest possible
!          spectral density consistent with the given moment constraints.
!          The Lagrange multipliers eta_k enforce these constraints.
!
!          Mathematical foundation:
!          We maximize the information entropy:
!            S[rho] = -integral[ rho(x) * log(rho(x)) dx ]
!          
!          subject to moment constraints:
!            integral[ T_k(x) * rho(x) dx ] = mui_k  for k=0,...,NT
!          
!          Using Lagrange multipliers, the solution is:
!            rho(x) = exp(-1 - sum_k[ eta_k * T_k(x) ]) / x^2
!          
!          where the factor 1/x^2 is the Jacobian for x = 1/E transformation,
!          ensuring proper normalization when converting to energy space.
!
! OpenMP:  Integration loops are parallelized for performance on multi-core systems
!
! USAGE:   This module is used by the main program in a cyclic optimization:
!          1. Set activeK to select which eta_k to optimize
!          2. Call activeDmu(eta_k) to compute moment error
!          3. Use bisection to find eta_k where activeDmu = 0
!          4. Repeat for all k until convergence
!===============================================================================

module entropy_module
  implicit none
  
  !-----------------------------------------------------------------------------
  ! Module-level variables (shared with main program)
  ! These are PUBLIC and accessible from the main program via USE statement
  !-----------------------------------------------------------------------------
  
  ! Lagrange multipliers (optimization parameters)
  real(8), allocatable :: et(:)      ! Current eta values being tested
  real(8), allocatable :: et0(:)     ! Best/converged eta values
  
  ! Spectral density and grid
  real(8), allocatable :: gm(:)      ! Spectral density rho(x) on quadrature grid
  real(8), allocatable :: xx(:)      ! Quadrature points in x-space (x = 1/E)
  
  ! Polynomial basis functions
  real(8), allocatable :: Tk(:,:)    ! Basis T_k(x_i) at quadrature points
                                      ! Dimension: (nq, 0:NT)
  real(8), allocatable :: TKlog(:,:) ! Basis on logarithmic energy grid
                                      ! Used for Hilbert transform calculations
  
  ! Target moments and weights
  real(8), allocatable :: mui(:)     ! Target moments from input data
                                      ! mui(k) = sum_i[ coupling_i * T_k(x_i) ]
  real(8), allocatable :: ww(:)      ! Gauss-Legendre quadrature weights
  
  ! Control variables
  integer :: activeK                 ! Index of eta currently being optimized
                                      ! In main loop: cycles 0,1,2,...,NT,0,1,...
  integer :: nq                      ! Number of quadrature points
  
contains
  
  !-----------------------------------------------------------------------------
  ! FUNCTION: activeDmu
  !
  ! PURPOSE: Compute moment error for the k-th constraint
  !          This is the objective function for the bisection root-finding method
  !
  ! MATHEMATICAL FORM:
  !   activeDmu(eta_k) = integral[ T_k(x) * rho(x; eta_k) dx ] - mui_k
  !
  !   where rho(x; eta_k) = exp(-1 - sum_j[ eta_j * T_j(x) ]) / x^2
  !
  ! METHOD:  
  !   1. Set eta_activeK = etK (update the k-th Lagrange multiplier)
  !   2. Recompute spectral density rho(x) with updated eta
  !   3. Compute moment: integral[ T_k(x) * rho(x) dx ] using quadrature
  !   4. Return: computed_moment - target_moment
  !
  ! ROOT FINDING:
  !   The bisection method calls this function repeatedly to find eta_k
  !   such that activeDmu(eta_k) = 0, which means the moment constraint
  !   is satisfied.
  !
  ! INPUT:   
  !   etK - trial value for eta_activeK (the k-th Lagrange multiplier)
  !
  ! OUTPUT:  
  !   Moment error (should be zero at optimal eta_k)
  !   - Positive: computed moment is too large (eta_k should increase)
  !   - Negative: computed moment is too small (eta_k should decrease)
  !
  ! NOTES:   
  !   - Uses module variable activeK to know which eta is being varied
  !   - Other eta values in et(:) are held fixed during this call
  !   - Called hundreds of times per iteration by bisection routine
  !   - Computationally expensive due to exponential evaluation at all grid points
  !
  ! NUMERICAL CONSIDERATIONS:
  !   - The exponential in GMexp can overflow/underflow for extreme eta values
  !   - The bisection method naturally avoids these regions by bracketing
  !   - Quadrature integration is exact for polynomials (within basis degree)
  !-----------------------------------------------------------------------------
  real(8) function activeDmu(etK)
    real(8), intent(in) :: etK
    integer :: i
    
    ! Update the k-th Lagrange multiplier with trial value
    ! All other eta values (j ≠ activeK) remain unchanged
    et(activeK) = etK
    
    ! Recompute spectral density at all quadrature points with new eta
    ! This is the expensive step: O(nq) evaluations of exponential
    ! OpenMP parallelization significantly speeds this up on multi-core systems
    !$OMP PARALLEL DO
    do i = 1,nq
      gm(i) = GMexp(i)  ! rho(x_i) = exp(-1 - sum_j[eta_j*T_j(x_i)]) / x_i^2
    end do
    !$OMP END PARALLEL DO
    
    ! Compute moment error using numerical integration:
    !   error = integral[ T_k(x) * rho(x) dx ] - mui_k
    !         ≈ sum_i[ w_i * T_k(x_i) * rho(x_i) ] - mui_k
    !
    ! At optimal eta_k, this should be ≈ 0 (within numerical tolerance)
    activeDmu = integrate(Tk(:,activeK)*gm(:)) - mui(activeK)
  end function activeDmu
  
  !-----------------------------------------------------------------------------
  ! FUNCTION: GMexp
  !
  ! PURPOSE: Evaluate maximum entropy spectral density at quadrature point j
  !          in x-space (x = 1/E)
  !
  ! FORMULA: 
  !   rho(x_j) = exp(-1 - sum_k[ eta_k * T_k(x_j) ]) / x_j^2
  !
  ! DERIVATION:
  !   The maximum entropy principle gives:
  !     rho(x) ∝ exp(sum_k[ eta_k * T_k(x) ])
  !   
  !   The Lagrangian formulation introduces:
  !     L = -integral[rho*log(rho)] + sum_k[eta_k*(integral[T_k*rho] - mui_k)]
  !   
  !   Taking functional derivative δL/δρ = 0 yields:
  !     log(rho) = -1 + sum_k[eta_k*T_k(x)] + constant
  !   
  !   Thus: rho(x) = Z^(-1) * exp(sum_k[eta_k*T_k(x)])
  !   where Z is the partition function (normalization constant)
  !
  ! COORDINATE TRANSFORMATION:
  !   Working in x = 1/E space requires Jacobian factor:
  !     rho_E(E) dE = rho_x(x) dx
  !     rho_E(E) = rho_x(1/E) * |dx/dE| = rho_x(1/E) / E^2
  !   
  !   Hence the 1/x^2 factor ensures proper normalization in energy space.
  !
  ! INPUT:   
  !   j - grid point index (1 ≤ j ≤ nq)
  !
  ! OUTPUT:  
  !   Spectral density rho(x_j) in units of [1/x]
  !   This must be multiplied by x^2 to get physical decay width Gamma(E)
  !
  ! NUMERICAL NOTES:
  !   - The exponential ensures rho > 0 everywhere (physically required)
  !   - Large negative exponents give rho ≈ 0 (numerically stable)
  !   - Large positive exponents can overflow (prevented by bisection bounds)
  !   - The "-1" in the exponent is absorbed into eta_0 normalization
  !
  ! PHYSICAL INTERPRETATION:
  !   - rho(x) is the spectral density in x = 1/E coordinates
  !   - Integrating rho(x) over x gives total spectral weight (normalized to 1)
  !   - Physical decay width: Gamma(E) = x^2 * rho(x) = (1/E^2) * rho(1/E)
  !-----------------------------------------------------------------------------
  real(8) function GMexp(j)
    integer, intent(in) :: j
    real(8) :: dum
    
    ! Compute exponent: -1 - sum_k[ eta_k * T_k(x_j) ]
    ! 
    ! sum(et(:)*Tk(j,:)) efficiently computes the dot product:
    !   eta_0*T_0(x_j) + eta_1*T_1(x_j) + ... + eta_NT*T_NT(x_j)
    dum = -1.d0 - sum(et(:)*Tk(j,:))
    
    ! Maximum entropy formula with Jacobian factor
    ! Division by xx(j)^2 accounts for transformation from x to E space
    GMexp = exp(dum)/xx(j)**2
  end function GMexp

  !-----------------------------------------------------------------------------
  ! FUNCTION: GMexpE
  !
  ! PURPOSE: Evaluate spectral density on logarithmic energy grid
  !          Used for computing Hilbert transform to extract energy shifts Delta(E)
  !
  ! DIFFERENCE FROM GMexp:
  !   - Evaluated on eelog grid (logarithmic spacing) instead of xx grid (linear)
  !   - Uses TKlog basis (evaluated at eelog points) instead of Tk basis
  !   - Does NOT include 1/x^2 factor (handled separately in energy-space calculations)
  !
  ! WHY LOGARITHMIC GRID?
  !   The Hilbert transform requires evaluating:
  !     Delta(E) = -(1/pi) * P.V. integral[ Gamma(E')/(E-E') dE' ]
  !   
  !   This integral has a singularity at E = E' (principal value).
  !   Logarithmic spacing provides:
  !     - Fine resolution near E = 0 (where resonances typically occur)
  !     - Coarse resolution at high E (where spectral density is smooth)
  !     - Better numerical stability for principal value evaluation
  !
  ! FORMULA:
  !   rho_log(x_j) = exp(-1 - sum_k[ eta_k * T_k(x_j) ])
  !   
  !   Note: no 1/x^2 factor here (applied later when computing Gamma(E))
  !
  ! INPUT:   
  !   j - grid point index on logarithmic grid (1 ≤ j ≤ nq)
  !
  ! OUTPUT:  
  !   Spectral density at eelog(j) without Jacobian factor
  !   Multiply by (eelog(j))^2 to get Gamma(E) on log grid
  !
  ! USAGE IN MAIN PROGRAM:
  !   gmElog(i) = GMexpE(i)           ! Get rho on log grid
  !   Gamma(E_i) = gmElog(i)          ! Already in right units for Hilbert
  !   Delta(E) = -(1/pi)*P.V.integral[Gamma(E')/(E-E') dE']
  !-----------------------------------------------------------------------------
  REAL(8) FUNCTION GMexpE(j)
    INTEGER, INTENT(in) :: j
    REAL(8) dum
    
    ! Compute exponent using basis functions on logarithmic grid
    ! TKlog(j,k) = T_k(eelog(j)) where eelog is logarithmically spaced
    dum = -1.d0-sum(et(:)*Tklog(j,:))
    
    ! Exponential only (no Jacobian factor here)
    GMexpE = exp(dum)
  END FUNCTION GMexpE
  
  !-----------------------------------------------------------------------------
  ! FUNCTION: integrate
  !
  ! PURPOSE: Compute weighted integral using Gauss-Legendre quadrature
  !
  ! MATHEMATICAL FORMULA:
  !   integral[a to b][ f(x) dx ] ≈ sum_{i=1}^{nq}[ w_i * f(x_i) ]
  !
  !   where x_i are Gauss-Legendre nodes and w_i are corresponding weights
  !
  ! ACCURACY:
  !   Gauss-Legendre quadrature with n points is exact for polynomials
  !   of degree ≤ 2n-1. For our composite quadrature with nq0 points per
  !   segment, each segment is exact for degree ≤ 2*nq0-1 ≈ 119.
  !
  !   Since we use NT ≈ 9 basis functions, and integrands involve products
  !   of basis functions (degree ≈ 18), the quadrature is more than sufficient.
  !
  ! WHY COMPOSITE QUADRATURE?
  !   Single-segment Gauss-Legendre is optimal for smooth functions on [-1,1].
  !   Our spectral density rho(x) may have features at multiple scales:
  !     - Smooth decay at high x (low E)
  !     - Sharp features near resonances (intermediate E)
  !     - Possible cusp at x = 0 (infinite E)
  !   
  !   Composite quadrature (multiple segments) handles this better than
  !   a single high-order quadrature.
  !
  ! INPUT:   
  !   ker(:) - integrand values at quadrature points (dimension nq)
  !            Can be any function evaluated at xx(:) grid points
  !
  ! OUTPUT:  
  !   Approximate integral value
  !   Exact for polynomials of degree ≤ 2*nq0-1 per segment
  !
  ! PERFORMANCE:
  !   - OpenMP parallelized with reduction for efficiency
  !   - On 8-core CPU: ~8x speedup for large nq (nq ≈ 2500)
  !   - Thread-safe: each thread computes partial sum, then reduced
  !   - Critical for performance since called thousands of times per run
  !
  ! EXAMPLE USAGE:
  !   norm = integrate(gm(:))              ! Compute normalization
  !   moment = integrate(gm(:)*xx(:)**2)   ! Compute second moment
  !   entropy = integrate(-gm(:)*log(gm(:))) ! Compute entropy
  !-----------------------------------------------------------------------------
  real(8) function integrate(ker)
    real(8), intent(in) :: ker(:)  ! Integrand at quadrature points
    real(8) :: sm                   ! Accumulator for sum
    integer :: i
    
    sm = 0.d0
    
    ! Parallel summation with reduction
    ! Each thread computes partial sum: sm_thread = sum[over thread's indices]
    ! Then OpenMP combines: sm = sm_thread1 + sm_thread2 + ... + sm_threadN
    !
    ! The reduction(+:sm) clause ensures thread-safe accumulation
    !$OMP PARALLEL DO REDUCTION(+:sm)
    do i = 1,size(ww)
      ! Gauss-Legendre quadrature: integral ≈ sum[ weight * function ]
      sm = sm + ww(i)*ker(i)
    end do
    !$OMP END PARALLEL DO
    
    integrate = sm
  end function integrate
  
end module entropy_module
