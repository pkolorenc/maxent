!===============================================================================
! PROGRAM: EntropyCode
!
! PURPOSE: Reconstructs energy-dependent decay width Gamma(E) from inverse 
!          spectral moments using Maximum Entropy Method (MEM).
!
! METHOD:  1. Evaluate spectral moments from discretized pseudocontinuum states
!          2. Maximize entropy subject to moment constraints using Lagrange
!             multipliers (eta parameters)
!          3. Iteratively solve for optimal eta values via bisection method
!          5. Reconstruct physical decay width Gamma(E) = x^2 * rho(x)
!             where x = 1/E and rho(x) is the spectral density on <0,1>
!
! INPUT:   coup.0001 - File containing:
!          - Line 1: working disctete state energy Ed (usually 0), true Ed wrt
!              neutral ground state
!          - Subsequent lines: index, energy of the pseudocont. state, coupling
!              to the discrete state
!
! OUTPUT:  - eta.last: Converged Lagrange multipliers (binary)
!          - File 2000: Polynomial expansion of Gamma(E)
!          - File 3000: Final entropy-maximizing Gamma(E)
!          - File 3001: level shift and Delta(E) and the width Gamma(E)
!          - File 9003: primitive spectral moments (input vs reconstructed)
!          - File 9006: working (Legendre) spectral moments (inp vs rec.)
!          - Screen: Convergence info, entropy, resonance parameters
!
! PHYSICS: The code works in x-space (x = 1/E) where rho(x) is normalized.
!          The physical decay width is Gamma(E) = (1/E^2) * rho(1/E).
!          At resonance energy Er, extracts Gamma(Er) (decay width) and 
!            Delta(Er) (energy shift) and via Hilbert transform of Gamma(E).
!
!===============================================================================

PROGRAM EntropyCode

  USE numerics_module
  USE entropy_module
  USE OMP_LIB

  IMPLICIT NONE

  ! Number of positive/logarithmic spectral moments to include (typically 1)
  INTEGER,PARAMETER :: Nlog = 1

  !-----------------------------------------------------------------------------
  ! Algorithm control parameters
  !-----------------------------------------------------------------------------
  INTEGER :: iterfact = 100000, &  ! Total iterations = iterfact * NT
             NT = 9, &             ! Number of polynomial basis functions T_k(x)
             nq0 = 60              ! Gauss-Legendre quadrature points per segment
  
  !-----------------------------------------------------------------------------
  ! Method selection flags
  !-----------------------------------------------------------------------------
  LOGICAL :: log_moment = .false., & ! Use logarithmic moment: mu_log = int[log(x)*rho(x)]
             inv_moment = .true., &  ! Use inverse moment: mu_{-1} = int[rho(x)/x] (default)
             ortho = .false., &      ! Apply Gram-Schmidt orthogonalization to basis
             full_X = .true., &      ! Extend integration domain to x=0 (infinite energy)
             devel = .false.         ! Enable developer mode (extra options)
  
  !-----------------------------------------------------------------------------
  ! Numerical parameters
  !-----------------------------------------------------------------------------
  REAL(8) :: tolerance = 1.d-09, &   ! Convergence criterion: sum(dmu^2) < tolerance
             E0shift = 0.d0, &       ! Manual shift to reference energy E0
             iter_scale = 1.00d0     ! Relaxation factor for eta updates (0 < alpha <= 1)

  CHARACTER(13) :: inpfile = "coup.0001"  ! Input coupling data filename
  CHARACTER(100) line                      ! Buffer for reading file header

  !-----------------------------------------------------------------------------
  ! Work variables
  !-----------------------------------------------------------------------------
  INTEGER ni,i,j,iou,akLast,itercheck, &
    itercount,ndum
  INTEGER,ALLOCATABLE :: count(:)          ! Failure counter for each eta_k
  
  ! Scalar work variables
  REAL(8) dx,dx0,e0,over,et1,dmu1,normscale,e_min,e_max, &
    dmu0_min,entropy_eei,low_xx_boundary, &
    high_xx_boundary,xx_range,renorm,dej,delta0,gamma0, &
    Ed,Eres,delres,gammares
  
  !-----------------------------------------------------------------------------
  ! Main arrays
  !-----------------------------------------------------------------------------
  REAL(8),ALLOCATABLE :: chWe(:), &        ! Chebyshev weights for input data points
    ee(:),Tke(:,:), &                      ! Energy grid and basis evaluated at data points
    dmu0(:),eta_old(:), &                  ! Moment errors and previous eta from file
    ww0(:),mu_inverse(:),et_best(:), &     ! Original quadrature weights, direct moments, best eta
    gmE(:),del(:),delker(:),eelog(:),gmElog(:)  ! E-space quantities for Hilbert transform
  
  LOGICAL dumlog,addLog  ! Flags for logic control

  !-----------------------------------------------------------------------------
  ! Namelist for runtime input parameters
  !-----------------------------------------------------------------------------
  namelist /entropy/ &
    NT, nq, nq0, iterfact, log_moment, inv_moment, tolerance, &
    legendre, inpfile, devel

  namelist /entdevel/ &
    chebyshev, ortho, E0shift, iter_scale, full_X

  !-----------------------------------------------------------------------------
  ! INITIALIZATION: Read runtime parameters from standard input
  !-----------------------------------------------------------------------------
  nq = 2500  ! Default number of quadrature points (will be adjusted)

  read(5,nml=entropy)
  write(6,nml=entropy)
  write(6,*) "==================="
  if (devel) then
    read(5,nml=entdevel)
    write(6,nml=entdevel)
    write(6,*) "==================="
  end if
  
  !-----------------------------------------------------------------------------
  ! Determine which additional moment constraint to use
  ! Options: inverse moment (1/x) or logarithmic moment log(x)
  ! These help constrain the spectral density at boundaries
  !-----------------------------------------------------------------------------
  if (inv_moment) log_moment = .false.  ! Inverse takes precedence
  if (log_moment.or.inv_moment) then
    addLog = .true.
  else
    addLog = .false.
  end if

  ! Report which entropy maximization variant is being used
  if (addLog) then
    if (inv_moment) then
      write(6,*) " trueEspace entropy max (first moment)"
    else
      write(6,*) " trueEspace entropy max (log)"
    end if
  else
    write(6,*) " trueEspace entropy max (0)"
  end if
  write(6,*) "==================="
  write(6,*)

  !-----------------------------------------------------------------------------
  ! STEP 1: Read input spectral moment data from file
  !
  ! File format:
  !   # E0 Ed                          (header: reference energy, discrete state energy)
  !   1  energy_shift_1  coupling_1   (data: index, energy, coupling strength)
  !   2  energy_shift_2  coupling_2
  !   ...
  !
  ! The coupling strengths represent the spectral weights from a discretized
  ! continuum calculation (e.g., from Feshbach projection or R-matrix theory)
  !-----------------------------------------------------------------------------
  write(6,*) "Reading input file:"
  call get_unit(iou)
  open(iou,file=inpfile,status="old")
  
  ! Read header line containing reference energy and discrete state energy
  read(iou,'(A)') line
  read(line(3:),*) e0,Ed
  Ed = Ed*Ha  ! Convert discrete state energy from Hartree to eV
  print*,e0,Ed
  
  ! Count number of data points in file
  i = 0
  ni = 0
  do while (i.eq.0)
    read(iou,*,iostat=i) j,dx,dx0
    if (i.eq.0) ni = ni+1
  end do
  
  allocate(ee(ni),gm(ni))
  rewind(iou)
  read(iou,*)  ! Skip header
  
  ! Read data points and convert to x = 1/E representation
  ! This transformation is key: working in x-space makes the problem better conditioned
  i = 0
  ni = 0
  e_min = 1.d8   ! Track energy range of input data
  e_max = -1.d8
  do while (i.eq.0)
    ni = ni+1
    read(iou,*,iostat=i) j,dx,dx0
    if (ni.eq.1) then
      ! First data point defines the reference energy
      e0 = 1.0d0-dx+E0shift
      print*,"e0",e0
    end if
    if (i.eq.0) then
      dx = dx+e0        ! Convert relative energy to absolute energy
      if (dx.lt.e_min) e_min = dx
      if (dx.gt.e_max) e_max = dx
      ee(ni) = 1.d0/dx  ! Transform to x-space: x = 1/E
      gm(ni) = dx0      ! Store coupling strength (will be squared later)
    end if
  end do
  close(iou)
  ni = size(gm)
  write(6,*) "... number of data points:",ni
  print*,ee(1),gm(1)    ! First data point (highest energy, smallest x)
  print*,ee(ni),gm(ni)  ! Last data point (lowest energy, largest x)

  ! Scale tolerance by number of basis functions
  tolerance = NT*tolerance

  !-----------------------------------------------------------------------------
  ! STEP 2: Setup integration domain in x-space [low_xx, high_xx]
  !
  ! The spectral density rho(x) will be reconstructed on this domain using
  ! composite Gauss-Legendre quadrature:
  !   - Divide range into segments
  !   - Use nq0 Gauss-Legendre points per segment
  !   - This handles multi-scale features better than single-segment quadrature
  !-----------------------------------------------------------------------------
  low_e_boundary = 1.d0+E0shift
  print*,"low_e_boundary:",low_e_boundary
  
  ! Set x-space boundaries
  low_xx_boundary = ee(ni)              ! Minimum x (corresponds to maximum energy)
  high_xx_boundary = 1.d0/low_e_boundary ! Maximum x (corresponds to minimum energy)
  
  ! Optional: extend integration to x=0 (infinite energy limit)
  ! This helps if spectral density has support at very high energies
  if (full_X) then
    low_xx_boundary = 0.0d0
  end if
  
  print*,"xx_boundary:",low_xx_boundary,high_xx_boundary
  xx_range = high_xx_boundary-low_xx_boundary
  print*,"xx_range:",xx_range

  ! Divide x-range into segments for composite quadrature
  ndum = int(nq/nq0)+1      ! Number of segments
  dx = xx_range/ndum        ! Segment width
  nq = ndum*nq0             ! Adjust total number of points
  print*,nq,ndum,nq0
  
  allocate(xx(nq),ww(nq),ww0(nq))
  
  ! Generate Gauss-Legendre points and weights for each segment
  ! Each segment [x_i, x_{i+1}] gets nq0 quadrature points
  do i = 0,ndum-1
    call GaussLeg(low_xx_boundary+i*dx,low_xx_boundary+(i+1.d0)*dx, &
      xx(i*nq0+1:(i+1)*nq0),ww0(i*nq0+1:(i+1)*nq0))
  end do
  print*,nq,xx(1),xx(nq)      ! Print x-range of quadrature
  print*,sum(ww0(:))          ! Should equal xx_range

  ww(:) = ww0(:)  ! Save original weights

  !-----------------------------------------------------------------------------
  ! STEP 3: Construct polynomial basis T_k(x) on quadrature grid
  !
  ! Options:
  !   - Legendre polynomials (default): orthogonal on [-1,1] with uniform weight
  !   - Chebyshev polynomials: orthogonal on [-1,1] with weight 1/sqrt(1-x^2)
  !   - Monomials x^k: simple but poorly conditioned
  !
  ! The choice affects numerical stability but not the physics (any complete
  ! basis can represent the spectral density)
  !-----------------------------------------------------------------------------
  if (addLog) then
    allocate(Tk(nq,0:NT+1))  ! Extra column for 1/x or log(x) basis function
  else
    allocate(Tk(nq,0:NT))
  end if
  
  ! Evaluate basis functions at all quadrature points
  do i = 1,nq
    call Trec01(xx(i),Tk(i,0:NT))  ! Standard polynomial basis
    if (addLog) then
      if (inv_moment) then
        Tk(i,NT+1) = 1.d0/xx(i)    ! Add 1/x to capture inverse moment
      else
        Tk(i,NT+1) = log(xx(i))    ! Add log(x) for logarithmic moment
      end if
    end if
  end do

  !-----------------------------------------------------------------------------
  ! STEP 4: Evaluate basis functions at input data points ee(:)
  !
  ! This is needed to compute the target moments:
  !   mui_k = sum_i [ coupling_i * T_k(x_i) ]
  !
  ! These moments are what the MEM solution must reproduce
  !-----------------------------------------------------------------------------
  allocate(chWe(ni))
  chWe(:) = 1.d0  ! Weight factors (could be Chebyshev weights if using that basis)
  
  if (addLog) then
    allocate(Tke(ni,0:NT+1))
  else
    allocate(Tke(ni,0:NT))
  end if
  
  ! Evaluate basis at input data points
  do i = 1,ni
    call Trec01(ee(i),Tke(i,0:NT))
    if (addLog) then
      if (inv_moment) then
        Tke(i,NT+1) = 1.d0/ee(i)
      else
        Tke(i,NT+1) = log(ee(i))
      end if
    end if
  end do

  !-----------------------------------------------------------------------------
  ! Setup logarithmic energy grid for Hilbert transform calculations
  !
  ! The Hilbert transform requires evaluating integrals with kernels like
  ! 1/(E-E'). Using a logarithmic grid provides better resolution where needed:
  !   - Fine spacing near E=0 (where resonances typically occur)
  !   - Coarse spacing at high energies (where spectral density is smooth)
  !-----------------------------------------------------------------------------
  allocate(eelog(nq))
  eelog(1) = 1.d0/xx(nq)    ! Highest energy point
  eelog(nq) = 1.d0/(xx(1)-xx(nq)+1.d0)  ! Lowest energy point
  
  ! Create logarithmic spacing
  dmu1 = log(eelog(1))
  dej = (log(eelog(nq))-dmu1)/(nq-1.d0)
  do i = 2,nq-1
    eelog(i) = exp(dmu1+(i-1)*dej)
  end do
  
  ! Convert back to x-space for consistency
  do i = 1,nq
    eelog(i) = 1.d0/eelog(i)
  end do
  
  ! Evaluate basis functions on logarithmic grid
  if (addlog) then
    allocate(TKlog(nq,0:NT+1))
  else
    allocate(TKlog(nq,0:NT))
  end if
  
  do i = 1,nq
    call Trec01(eelog(i),Tklog(i,0:NT))
    if (addLog) then
      if (inv_moment) then
        Tklog(i,NT+1) = 1.d0/eelog(i)
      else
        Tklog(i,NT+1) = log(eelog(i))
      end if
    end if
  end do

  !-----------------------------------------------------------------------------
  ! STEP 5: Compute spectral moments from input data
  !
  ! Two types of moments:
  !   1. Direct moments: mu_j = sum_i [ gm(i) * x_i^j ]
  !      These are used for verification (compare input vs reconstructed)
  !
  !   2. Projected moments: mui_k = sum_i [ gm(i) * T_k(x_i) ]
  !      These are the actual constraints for MEM
  !
  ! Physical interpretation:
  !   - mu_0: total spectral weight (normalization)
  !   - mu_1: first moment (related to centroid energy)
  !   - Higher moments: shape information
  !-----------------------------------------------------------------------------
  allocate(mu_inverse(-2:50))
  mu_inverse = 0.d0
  
  do i = 1,ni
    ! Convert coupling strength to physical units:
    ! Factor breakdown:
    !   - Ha*1000: Hartree to meV conversion
    !   - 2*pi: frequency to energy conversion
    !   - ^2: |coupling|^2 gives transition rate
    gm(i) = Ha*1000*2.d0*pi*gm(i)**2
    
    ! Compute direct moments mu_j = sum[ gm(i) * x_i^j ] for j = -2 to 50
    do j = -2,50
      mu_inverse(j) = mu_inverse(j)+gm(i)*ee(i)**j
    end do
  end do

  !-----------------------------------------------------------------------------
  ! STEP 6: Optional Gram-Schmidt orthogonalization of basis
  !
  ! This creates an orthonormal basis specifically for the given quadrature:
  !   <T_i, T_j> = integral[ T_i(x) * T_j(x) dx ] = delta_ij
  !
  ! Advantages:
  !   - Better numerical conditioning
  !   - Diagonal moment matrix
  !
  ! Disadvantages:
  !   - Loses analytic properties of Legendre/Chebyshev polynomials
  !   - More computational cost
  !-----------------------------------------------------------------------------
  if (ortho) then
    do i = 0,NT
      ! Orthogonalize T_i against all previous T_j (j < i)
      do j = 0,i-1
        over = integrate(Tk(:,i)*Tk(:,j))  ! Compute overlap
        Tk(:,i) = Tk(:,i)-over*Tk(:,j)     ! Subtract projection
        Tke(:,i) = Tke(:,i)-over*Tke(:,j)  ! Also at data points
        Tklog(:,i) = Tklog(:,i)-over*Tklog(:,j)  ! And log grid
      end do
      ! Normalize T_i
      over = integrate(Tk(:,i)**2)
      Tk(:,i) = Tk(:,i)/sqrt(over)
      Tke(:,i) = Tke(:,i)/sqrt(over)
      Tklog(:,i) = Tklog(:,i)/sqrt(over)
    end do
  end if

  !-----------------------------------------------------------------------------
  ! STEP 7: Normalize each polynomial basis function
  !
  ! This is done regardless of orthogonalization to ensure:
  !   integral[ T_k(x)^2 dx ] = 1
  !
  ! Normalization improves numerical stability and makes eta values
  ! have comparable magnitudes
  !-----------------------------------------------------------------------------
  do i = 0,NT
    over = integrate(Tk(:,i)**2)
    Tk(:,i) = Tk(:,i)/sqrt(over)
    Tke(:,i) = Tke(:,i)/sqrt(over)
    Tklog(:,i) = Tklog(:,i)/sqrt(over)
  end do

  ! Adjust NT to include logarithmic/inverse moment if used
  if (addLog) NT = NT+Nlog

  !-----------------------------------------------------------------------------
  ! STEP 8: Project input data onto polynomial basis to get target moments
  !
  ! Compute: mui(k) = sum_i[ gm(i) * T_k(x_i) ]
  !
  ! These are the constraints the MEM solution must satisfy:
  !   integral[ T_k(x) * rho_MEM(x) dx ] = mui(k)  for k=0,...,NT
  !
  ! The maximum entropy principle then gives the unique rho(x) with maximum
  ! entropy that satisfies these constraints
  !-----------------------------------------------------------------------------
  allocate(mui(0:NT))
  mui(:) = 0.d0
  
  do i = 0,NT
    ! Compute moment as weighted sum over data points
    do j = size(gm),1,-1
      mui(i) = mui(i)+chWe(j)*gm(j)*Tke(j,i)
    end do
    
    ! Ensure positivity by flipping sign if needed
    ! (This is a convention choice; could also flip the basis function)
    if (mui(i).lt.0.d0) then
      mui(i) = -mui(i)
      Tke(:,i) = -Tke(:,i)
      Tk(:,i) = -Tk(:,i)
      Tklog(:,i) = -Tklog(:,i)
    end if
  end do

  write(6,*)
  print*,"mui raw =",mui
  
  ! Save normalization constant and normalize moments to unit total weight
  normscale = mui(0)
  write(6,*)
  print*,"scale = ",normscale
  mui(:) = mui(:)/normscale
  
  deallocate(ee,gm,chWe,Tke)
  write(6,*)
  print*,"mui scaled =",mui
  write(6,*)

  !-----------------------------------------------------------------------------
  ! Evaluate direct polynomial expansion (before MEM optimization)
  !
  ! This is simply: rho_poly(x) = sum_k [ mui_k * T_k(x) ]
  !
  ! This satisfies the moment constraints by construction, but may have
  ! unphysical oscillations. The MEM optimization will smooth it out while
  ! maintaining the moment constraints.
  !-----------------------------------------------------------------------------
  dx = 1000.0d0
  dx0 = 0.d0
  allocate(ee(nq),gm(nq))
  
  if (addLog) then
    j = NT-Nlog  ! Exclude log/inverse moment from polynomial expansion
  else
    j = NT
  end if
  
  ! Compute polynomial expansion at each quadrature point
  do i = 1,nq
    eelog(i) = 1.d0/eelog(i)-e0  ! Convert to energy relative to E0
    ee(i) = 1.d0/xx(i)-e0
    gm(i) = sum(mui(:j)*Tk(i,:j))  ! rho_poly(x) = sum[ mui_k * T_k(x) ]
    if (gm(i).lt.0.d0) gm(i) = 0.d0  ! Enforce positivity
  end do
  
  ! Verify and enforce normalization
  dx = sum(ww0(:)*gm(:))
  write(6,*) "Gm_poly norm: ",dx
  gm(:) = gm(:)/dx  ! Renormalize
  dx = sum(ww0(:)*gm(:))
  write(6,*) "Gm_poly renorm: ",dx
  
  ! Output polynomial expansion to file 2000
  ! Columns: E, Gamma(E), x, rho(x)
  do i = nq,1,-1
    write(2000,'(4ES20.10)') ee(i),normscale*gm(i)*xx(i)**2,xx(i),gm(i)
    ! Find Gamma(E=0) for reference
    if (abs(ee(i)).lt.dx) then
      dx = abs(ee(i))
      dx0 = gm(i)*xx(i)**2
    end if
  end do

  write(6,*) "Gm(0)_poly = ",dx0*normscale
  
  ! Compute entropy of polynomial expansion: S = -integral[ rho * log(rho) ]
  ee(:) = 0.d0
  do i = 1,nq
    if (gm(i).gt.0.d0) ee(i) = -gm(i)*log(gm(i))
  end do
  write(6,*) "entropy_poly = ",sum(ww0(:)*ee(:))

  !-----------------------------------------------------------------------------
  ! STEP 9: Initialize Lagrange multipliers eta_k
  !
  ! The MEM solution has the form:
  !   rho(x) = exp(-1 - sum_k[ eta_k * T_k(x) ]) / x^2
  !
  ! The eta_k are Lagrange multipliers that enforce the moment constraints.
  ! Finding optimal eta is the core computational task.
  !
  ! Initialization strategy:
  !   1. Try to load eta from previous run (warm start)
  !   2. If not available, use random seed near zero
  !   3. Special handling for eta_0 (normalization) and eta_NT (log/inverse moment)
  !-----------------------------------------------------------------------------
  allocate(et0(0:NT),et(0:NT),dmu0(0:NT))
  
  ! Generate random initial guess
  call random_seed()
  call random_number(et0(:))
  do i = 0,NT-1
    et0(i) = 2.d0*et0(i)-1.d0  ! Map [0,1] to [-1,1]
    et0(i) = 1.e-2*et0(i)      ! Scale to small values
  end do
  
  ! Special initialization for key parameters
  et0(0) = 1.d0  ! eta_0 controls normalization (start near 1)
  if (addLog) then
    if (inv_moment) then
      et0(NT) = 0.1d0   ! Inverse moment eta
    else
      et0(NT) = 1.1d0   ! Log moment eta
    end if
  end if

  ! Try loading previous eta values from binary file
  call get_unit(iou)
  write(6,*)
  open(iou,file="eta.last",form="unformatted",access="sequential", &
    status="old",iostat=i)
  if (i.ne.0) then
    write(6,*) "... starting from random seed"
  else
    print*,"previous eta loaded"
    read(iou) activeK,dumlog
    print*,activeK,dumlog
    allocate(eta_old(0:activeK))
    read(iou) eta_old(:)
    print*,"eta_old:"
    print*,eta_old
    et0(:) = 0.d0

    ! Match old eta to new basis size (handles case where NT changed)
    if (addLog) then
      if (dumlog) then
        ! Both old and new use log/inverse moment
        do i = 0,min(NT-Nlog,activeK-Nlog)
          et0(i) = eta_old(i)
        end do
        et0(NT) = eta_old(activeK)
        if (Nlog.eq.2) et0(NT-1) = eta_old(activeK-1)
      else
        ! Old didn't use log/inverse, new does
        do i = 0,min(NT-Nlog,activeK)
          et0(i) = eta_old(i)
        end do
      end if
    else
      if (dumlog) then
        ! Old used log/inverse, new doesn't
        do i = 0,min(NT,activeK-Nlog)
          et0(i) = eta_old(i)
        end do
      else
        ! Neither uses log/inverse
        do i = 0,min(NT,activeK)
          et0(i) = eta_old(i)
        end do
      end if
    end if
    write(6,*) "... starting from last eta"
  end if
  print*,et0

  !-----------------------------------------------------------------------------
  ! STEP 10: Evaluate initial moment errors with eta_0
  !
  ! For each constraint k, compute:
  !   dmu_k = integral[ T_k(x) * rho(x;eta_0) dx ] - mui_k
  !
  ! The goal of optimization is to make all dmu_k ≈ 0
  ! We minimize the objective function: sum_k[ dmu_k^2 ]
  !-----------------------------------------------------------------------------
  dmu0 = 0.d0
  et(:) = et0(:)
  do activeK = 0,NT
    dmu0(activeK) = activeDmu(et0(activeK))
    print*,activeK,dmu0(activeK)**2
  end do
  dmu0_min = sum(dmu0(:)**2)  ! Initial objective function value
  write(6,*)
  print*,"dmu2 =",dmu0_min

  !-----------------------------------------------------------------------------
  ! STEP 11: Main iteration loop to optimize eta parameters
  !
  ! Algorithm: Cyclic Coordinate Descent with Bisection
  !
  ! For each iteration:
  !   1. Select one eta_k to optimize (cycle through k = 0,1,...,NT)
  !   2. Fix all other eta_j (j ≠ k)
  !   3. Find eta_k that makes dmu_k(eta_k) = 0 using bisection
  !   4. Update eta_k with relaxation
  !   5. Periodically check global convergence
  !
  ! Convergence criteria:
  !   - Primary: sum(dmu^2) < tolerance
  !   - Secondary: No improvement for 20*NT checks (stalled)
  !   - Failure: Individual eta not converging after NT+2 attempts
  !
  ! This approach is robust but can be slow. The key insight is that
  ! optimizing one eta at a time is numerically stable, even though
  ! the full optimization problem is nonlinear and high-dimensional.
  !-----------------------------------------------------------------------------
  allocate(count(0:NT))
  count = 0
  write(6,*) "Iterations"
  aKlast = -1
  allocate(et_best(0:size(et0)-1))
  et_best = et0
  itercheck = 0
  itercount = 0
  
  ! Main iteration loop: cycle through all eta parameters many times
  do j = iterfact*NT-1,0,-1
    activeK = mod(j,NT+1)  ! Cycle: 0,1,2,...,NT,0,1,2,...,NT,...
    
    ! Set search range for bisection method
    ! The range adapts to current eta value to avoid overshooting
    if (abs(et0(activeK)).gt.0.d0) then
      dx = abs(100.d0*et0(activeK))  ! Search within ±100*|eta|
    else
      dx = 100.0  ! Default range for eta ≈ 0
    end if
    
    ! Find eta_k that zeros activeDmu(eta_k) using bisection
    ! This solves: integral[ T_k(x) * rho(x;eta_k) dx ] = mui_k
    call bisection(activeDmu,et0(activeK)-dx,et0(activeK)+dx,1.d-10,50,et1,dumlog)
    
    if (dumlog) then
      ! Bisection succeeded: update eta with relaxation
      ! Relaxation (iter_scale < 1) prevents oscillations
      et0(activeK) = et0(activeK) + iter_scale*(et1-et0(activeK))
      count(activeK) = 0
    else
      ! Bisection failed: count failure and potentially exit
      count(activeK) = count(activeK)+1
      print*,j,i
      if (count(activeK).gt.NT+2) then
        write(6,*) "not converging"
        exit
      end if
    end if
    
    !---------------------------------------------------------------------------
    ! Periodic convergence check (every 100*NT+1 iterations)
    !---------------------------------------------------------------------------
    if (mod(j,100*NT+1).eq.itercheck) then
      itercheck = mod(itercheck+1,NT+1)
      aKlast = activeK
      
      ! Recompute all moment errors with current eta
      do activeK = 0,NT
        dmu0(activeK) = activeDmu(et0(activeK))
      end do
      dx = sum(dmu0(:)**2)  ! Current objective function value
      write(6,'(I9," ",I2," ",ES10.3," ",ES10.3," ",I4)') &
        j,aKlast,dx,dx/tolerance,itercount
      
      ! Track best eta found so far
      if (dx.lt.dmu0_min) then
        ! Significant improvement: reset stall counter
        if ((dmu0_min-dx)/dx.gt.1.d-5) itercount = 0
        dmu0_min = dx
        et_best = et0
      else
        ! No improvement: increment stall counter
        itercount = itercount+1
      end if
      
      ! Check convergence criteria
      if (dx.lt.tolerance) then
        write(6,*)
        write(6,*) " +++ tolerance reached +++"
        exit
      end if
      if (itercount.ge.20*NT) then
        write(6,*)
        write(6,*) " +++ stalled +++"
        exit
      end if
    end if
  end do

  ! Use best eta found (in case last iteration wasn't the best)
  if (dmu0_min.lt.1.d15) et0 = et_best
  
  ! Final moment error evaluation
  do activeK = 0,NT
    dmu0(activeK) = activeDmu(et0(activeK))
  end do
  print*,"=========="
  print*,"dmu2 =",sum(dmu0(:)**2)
  write(6,*)
  print*,"et0 =",et0
  write(6,*)
  print*,"dmu =",dmu0

  !-----------------------------------------------------------------------------
  ! Output final results
  !-----------------------------------------------------------------------------
  write(6,*) "------------------------------------------"
  write(6,*)
  write(6,*) "Final result:"
  write(6,*) "============="

  ! Report final moment errors (absolute and relative)
  do activeK = 0,NT
    dmu0(activeK) = activeDmu(et0(activeK))
  end do
  write(6,*)
  print*,"dmu2 (final) =",sum(dmu0(:)**2)
  write(6,*)
  print*,"et0(normalized) ="
  print*,et0
  write(6,*)
  do activeK = 0,NT
    dmu0(activeK) = abs(dmu0(activeK)/mui(activeK))
  end do
  print*,"dmu (rel) =",dmu0
  write(6,*)

  !-----------------------------------------------------------------------------
  ! Compute final spectral density rho(x) from optimized eta
  !
  ! The MEM solution is:
  !   rho(x) = exp(-1 - sum_k[ eta_k * T_k(x) ]) / x^2
  !
  ! This form guarantees:
  !   1. Positivity: rho(x) > 0 everywhere (exponential)
  !   2. Moment constraints: integral[ T_k * rho ] = mui_k (via eta optimization)
  !   3. Maximum entropy: S = -integral[ rho * log(rho) ] is maximized
  !-----------------------------------------------------------------------------
  et = et0
  dx0 = 0.d0
  if (allocated(ee)) deallocate(ee,gm)
  allocate(ee(nq),gm(nq))
  
  do i = 1,nq
    ee(i) = 1.d0/xx(i)-e0  ! Convert x back to energy
    gm(i) = GMexp(i)       ! Evaluate MEM spectral density
  end do

  !---------------------------------------------------------------------------
  ! Verify moments are reproduced correctly
  ! Output to file 9003: comparison of input vs reconstructed moments
  !---------------------------------------------------------------------------
  do j = -2,50
    dx = normscale*integrate(gm(:)*xx(:)**j)
    write(9003,'(I3," ",3ES25.16)') j, &
      mu_inverse(j), &                ! Input moment
      dx, &                           ! Reconstructed moment
      abs((dx-mu_inverse(j))/mu_inverse(j))  ! Relative error
  end do

  !-----------------------------------------------------------------------------
  ! Compute physical decay width Gamma(E)
  !
  ! Relation: Gamma(E) = (1/E^2) * rho(1/E) = x^2 * rho(x)
  !
  ! Physical interpretation:
  !   - Gamma(E) is the energy-dependent decay rate of the resonance
  !   - Units: meV (milli-electron volts)
  !   - Integrated over E: gives total decay rate
  !-----------------------------------------------------------------------------
  dx = 1000.0d0
  dumlog = .true.
  allocate(gmE(size(ee)))
  
  do i = nq,1,-1
    gmE(i) = normscale*gm(i)*xx(i)**2  ! Gamma(E) = x^2 * rho(x)
    write(3000,'(4ES25.10E3)') ee(i),gmE(i),xx(i),normscale*gm(i)
    
    ! Find Gamma(E=0) by linear extrapolation
    if (dumlog) then
      if (ee(i).gt.0.d0) then
        dumlog = .false.
        dx = gmE(i)-gmE(i+1)
        dx0 = gmE(i+1)-dx*ee(i+1)/(ee(i)-ee(i+1))
      end if
    end if
  end do

  !-----------------------------------------------------------------------------
  ! Compute energy-dependent level shift Delta(E) via Hilbert transform
  !
  ! The Kramers-Kronig relation connects real and imaginary parts of
  ! the self-energy. For resonances:
  !
  !   Delta(E) = -(1/pi) * P.V. integral[ Gamma(E')/(E-E') dE' ]
  !
  ! where P.V. denotes principal value (exclude singularity at E=E')
  !
  ! Physical interpretation:
  !   - Delta(E) is the energy shift due to coupling to continuum
  !   - Renormalizes the discrete state energy: E_res = E_d + Delta(E_res)
  !   - Units: meV
  !-----------------------------------------------------------------------------
  ni = size(eelog)
  allocate(del(ni),delker(ni),gmElog(ni))
  
  ! Evaluate Gamma(E) on logarithmic grid for better resolution
  do i = 1,ni
    gmElog(i) = GMexpE(i)
  end do
  gmElog(:) = normscale*gmElog(:)
  
  del = 0.d0
  
  ! For each energy point E_j, compute Hilbert transform
  do j = 1+2,ni-2  ! Skip endpoints to avoid boundary issues
    delker(:) = 0.d0
    
    do i = 1,ni
      if (i.ne.j) then
        ! Regular points: kernel = Gamma(E_i)/(E_j - E_i)
        delker(i) = gmElog(i)/(eelog(j)-eelog(i))
      else
        ! Singular point: use finite difference approximation
        ! This implements the principal value prescription
        delker(j) = gmElog(j+1)-gmElog(j-1)
      end if
    end do
    del(j) = delker(j)
    
    ! Integrate from left side (E < E_j)
    dmu1 = 0.d0
    do i = 1,j-2
      dmu1 = dmu1+0.5d0*(delker(i)+delker(i+1))*(eelog(i+1)-eelog(i))
    end do
    del(j) = del(j)+dmu1
    
    ! Integrate from right side (E > E_j)
    dmu1 = 0.d0
    do i = ni,j+2,-1
      dmu1 = dmu1+0.5d0*(delker(i-1)+delker(i))*(eelog(i)-eelog(i-1))
    end do
    del(j) = del(j)+dmu1
  end do
  
  ! Apply prefactor: Delta(E) = -(1/pi) * integral
  del(:) = -del(:)*0.5d0/pi
  deallocate(delker)
  
  !---------------------------------------------------------------------------
  ! Extract Delta(Ed) and Gamma(Ed) at discrete state energy
  ! Use linear interpolation to find values at Ed
  !---------------------------------------------------------------------------
  dumlog = .true.
  do j = 1,ni
    write(3001,'(3ES25.10E3)') eelog(j),gmElog(j),del(j)
    if (dumlog) then
      if (eelog(j).gt.0.d0) then
        dumlog = .false.
        ! Linear interpolation to E=0
        delta0 = del(j-1)-(del(j)-del(j-1))*eelog(j-1)/(eelog(j)-eelog(j-1))
        gamma0 = gmElog(j-1)-(gmElog(j)-gmElog(j-1))*eelog(j-1)/ &
          (eelog(j)-eelog(j-1))
      end if
    end if
  end do

  !-----------------------------------------------------------------------------
  ! Find resonance energy Er where Delta(Er) = Er (self-consistency)
  !
  ! The physical resonance occurs where the self-consistent equation is satisfied:
  !   E_r = E_d + Delta(E_r)
  !
  ! Rearranging: E_r - Delta(E_r) = 0
  !
  ! We solve this by finding where the curve E vs Delta(E)/1000 crosses E=Delta
  !-----------------------------------------------------------------------------
  eelog(:) = eelog(:)*Ha  ! Convert to Hartree
  dumlog = .true.
  
  do j = 1,ni
    if (dumlog) then
      if (eelog(j).gt.del(j)/1000.d0) then
        dumlog = .false.
        ! Linear interpolation to find crossing point
        Eres = del(j-1)/1000.d0*(eelog(j)-eelog(j-1))/ &
          ((del(j-1)-del(j))/1000.d0+eelog(j)-eelog(j-1))
        delres = del(j-1)+(del(j)-del(j-1))/(eelog(j)-eelog(j-1))* &
          (Eres-eelog(j-1))
        gammares = gmElog(j-1)+(gmElog(j)-gmElog(j-1))/(eelog(j)-eelog(j-1))* &
          (Eres-eelog(j-1))
      end if
    end if
  end do
  Eres = Eres+Ed  ! Add back discrete state energy reference
  
  deallocate(del,eelog,gmElog)

  !-----------------------------------------------------------------------------
  ! Output moment reproduction quality to file 9006
  ! Shows how well the MEM solution satisfies the moment constraints
  !-----------------------------------------------------------------------------
  if (addLog) then
    activeK = NT
    dmu0(activeK) = activeDmu(et0(activeK))
    write(9006,'(I3," ",4ES25.16)') -1,mui(activeK)*normscale, &
      normscale*(dmu0(activeK)+mui(activeK)),abs(dmu0(activeK)/mui(activeK)), &
      mui(activeK)
    i = NT-1
  else
    i = NT
  end if
  
  do activeK = 0,i
    dmu0(activeK) = activeDmu(et0(activeK))
    write(9006,'(I3," ",4ES25.16)') activeK,mui(activeK)*normscale, &
      normscale*(dmu0(activeK)+mui(activeK)),abs(dmu0(activeK)/mui(activeK)), &
      mui(activeK)
  end do
  write(6,*)
  print*,"dmu2 (final) =",sum(dmu0(:)**2)

  !---------------------------------------------------------------------------
  ! Final normalization check
  ! The MEM solution should integrate to 1 by construction
  !---------------------------------------------------------------------------
  dx = sum(ww0(:)*gm(:))
  write(6,*) "Gm_ent norm: ",abs(1.d0-dx)
  renorm = 1.d0/dx
  gm(:) = gm(:)*renorm
  dx = sum(ww0(:)*gm(:))
  write(6,*) "Gm_ent renorm: ",abs(1.d0-dx)

  write(6,*) "Gm(0)_ent = ",dx0

  !-----------------------------------------------------------------------------
  ! Print physical results to screen
  !-----------------------------------------------------------------------------
  write(6,*)
  write(6,*) "Ed [eV] = ",Ed
  write(6,*) "Delta(Ed) [meV] = ",delta0
  write(6,*) "Gamma(Ed) [meV] = ",gamma0

  write(6,*)
  write(6,*) "Er [eV] = ", Eres
  write(6,*) "Delta(Er) [meV] = ",delres
  write(6,*) "Gamma(Er) [meV] = ",gammares

  !-----------------------------------------------------------------------------
  ! Compute entropy of final solution
  !
  ! Two definitions:
  !   S_x = -integral[ rho(x) * log(rho(x)) dx ]      (in x-space)
  !   S_E = -integral[ rho(x) * log(x^2*rho(x)) dx ]  (in E-space)
  !
  ! The MEM maximizes S_x subject to constraints
  ! S_E accounts for Jacobian of transformation x = 1/E
  !-----------------------------------------------------------------------------
  ee(:) = 0.d0
  do i = 1,nq
    if (gm(i).gt.0.d0) ee(i) = -gm(i)*log(gm(i))
  end do
  entropy_eei = integrate(ee(:))  ! S_x
  
  do i = 1,nq
    if (gm(i).gt.0.d0) ee(i) = -gm(i)*log(xx(i)**2*gm(i))
  end do
  write(6,*)
  write(6,*) "entropy_ent = ",entropy_eei,integrate(ee(:))
  write(6,*) " ... S_x, S_e"

  !-----------------------------------------------------------------------------
  ! Save converged eta parameters for next run (warm start)
  ! This enables iterative refinement: run with increasing NT
  !-----------------------------------------------------------------------------
  call get_unit(iou)
  open(iou,file="eta.last",form="unformatted",access="sequential", &
    status="unknown")
  write(iou) NT,addLog
  write(iou) et0
  close(iou)

  ! Adjust eta_0 to account for normalization constant
  ! This makes eta.last contain the full solution including normalization
  et0(0) = et0(0)+log(normscale)/Tk(1,0)
  write(6,*)
  print*,"et0 = "
  print*,et0(:)
  write(6,*) "NB: eta.last contains exponents of 1-normalized probability density"

END PROGRAM EntropyCode
