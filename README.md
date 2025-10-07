# Maximum Entropy Method for the reconstruction of the decay width function from spectral moments

&nbsp;

[![DOI](https://zenodo.org/badge/1071499585.svg)](https://doi.org/10.5281/zenodo.17287708)

&nbsp;

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Output Files](#output-files)


---

## Overview

This code reconstructs energy-dependent decay width function Γ(E) from couplings to a discretized continuum using the **Maximum Entropy Method (MEM)** for an inverse moment problem.

&nbsp;

### What does it do?

Given a finite set of coupling strengths at discrete pseudo-continuum energies (from quantum chemistry or scattering calculations), the code:
- Evaluates spectral moments with respect to Legendre polynomials (or primitive power functions)
- Reconstructs the continuous spectral density ρ(x) in x = 1/E space
- Extracts the physical decay width Γ(E) = (1/E²)ρ(1/E)
- Identifies resonance parameters (position Er, decay width Γr)

&nbsp;


### The Maximum Entropy Principle

**Problem**: Given N spectral moments {μ₀, μ₁, ..., μₙ}, find the spectral density ρ(x).

**Solution**: Maximize information entropy S = -∫ ρ log(ρ) dx subject to moment constraints.

**Result**: Unique, smooth, positive spectral density:
```
ρ(x) = exp(-1 - Σₖ ηₖ Tₖ(x)) / x²
```

where ηₖ are Lagrange multipliers and Tₖ(x) are (orthogonal) polynomials.

&nbsp;

### Algorithm Flow


```
1. Read input data (coupling strengths at discrete energies)
   ↓
2. Transform to x-space: x = 1/E
   ↓
3. Compute target moments: μₖ = Σᵢ[coupling_i × Tₖ(xᵢ)]
   ↓
4. Initialize Lagrange multipliers ηₖ (random or from file)
   ↓
5. Main optimization loop:
   For iter = 1 to iterfact×NT:
     a) Select k = mod(iter, NT+1)
     b) Fix all ηⱼ (j≠k)
     c) Solve for ηₖ: ∫Tₖ(x)ρ(x;ηₖ)dx = μₖ using bisection
     d) Update ηₖ with relaxation
     e) Check convergence every 100×NT iterations
   ↓
6. Compute final spectral density: ρ(x) = exp(-1 - Σₖ ηₖTₖ(x))/x²
   ↓
7. Transform to physical observables:
   - Decay width: Γ(E) = x²ρ(x)
   - Energy shift: Δ(E) via Hilbert transform
   - Resonance: solve Er = Ed + Δ(Er)
   ↓
8. Output results and save ηₖ for next run
```

&nbsp;


---

## Installation
&nbsp;

### Requirements

- **Compiler**: Intel Fortran (ifx/ifort) or GNU Fortran (gfortran)
- **OpenMP**: For parallelization (standard with modern compilers)
- **No external libraries required**


### Build outputs
- `entropy.exe`: Main program
- `readEta.exe`: Utility to read binary eta.last file


### Running the Code
- complete pipeline is available in the `RUN` directory with expected output in  the `TEST.tar.gz` archive

&nbsp;

---

## Output Files
&nbsp;

### Binary file (restart)
| File | Description | Format |
|------|-------------|--------|
| `eta.last` | Converged Lagrange multipliers | Unformatted binary |

&nbsp;

### ASCII output files (fort.*)
| File | Content | Columns |
|------|---------|---------|
| `fort.2000` | Polynomial expansion Γ(E) | E, Γ(E), x, ρ(x) |
| `fort.3000` | MEM solution Γ(E) | E, Γ(E), x, ρ(x) |
| `fort.3001` | Energy shifts and widths | E, Γ(E), Δ(E) |
| `fort.9003` | Moment verification (primitive moments) | j, μⱼ(input), μⱼ(recon), error |
| `fort.9006` | Moment verification (working moments - Legendre/Cheb/prim) | k, μₖ, reconstructed, rel. error |

&nbsp;

---


## License


This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this code in published work, please cite:

```
A. Chalúpek and P. Kolorenč: Entropy solution to the moment problem in molecular decay width calculations. [to be published].
```

