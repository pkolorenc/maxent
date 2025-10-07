# MaxEnt Algorithm Pipeline

&nbsp;

## Overview

This pipeline executes the maximum entropy algorithm to solve the moment problem for the decay width. It runs calculations across a range of moment orders (M=0-20 by default), extracts key physical quantities, and computes averaged decay width functions functions with uncertainties.

&nbsp;

## Directory Structure

```

├── runMaxEnt.sh           # Main driver script
├── aux/
│   ├── collect.sh         # Results extraction
│   ├── inp_subs.py        # Input file parameter substitution
│   ├── get_max_dmu.py     # spectral moments convergence analysis
│   └── maxentAverage.py   # Gamma(E) function averaging
├── entropy.input          # Input parameters template
├── coup.0001              # input data
└── ../
    ├── entropy.exe        # Main MaxEnt executable
    └── readEta.exe        # Eta file reader
```
&nbsp;

## Quick Start

&nbsp;

### Prerequisites and input

&nbsp;

1. Compiled executables:
   - `entropy.exe` (main MaxEnt solver)
   - `readEta.exe` (eta file parser)

2. Input files in working directory:
   - `entropy.input` (input file template for the main solver)
        * NT = 00               ! M ... is replaced by the driver script for each run
        * NQ = 2000             ! size of the working quadrature on the <0,1> x-interval
        * ITERFACT = 40000      ! max number of iterations = NT*ITERFACT
        * TOLERANCE = 1.e-9     ! convergence criterion for D(lambda) is NT*TOLERANCE
        * LEGENDRE = T          ! compute spectral moments with respect to Legendre polynomials
        * INV_MOMENT = T        ! add x^{-1} spectral moment (DEFAULT)
        * LOG_MOMENT = F        ! add logarithmic moment (INV_MOMENT has priority if both T)
        * INPFILE = "coup.0001" ! name of the input data file
        
   - `coup.0001` (input data)
        * on the first line, it contains the working discrete state energy (typicall 0.0 if the pseudo-continuum energy is given relative to Ed) and the actual discrete state energy relative to neutral ground state (all in Ha)
        * the following lines contain coupling data in Ha: `i  eps_i-Ed   <d|H|chi_i>`
   - `eta.last` (optional - used to initialize the Lagrange multipliers)

3. Environment setup:
   - `./envir.sh` should configure the environment: Intel compiler, OpenMP threads, stack size, ...

   &nbsp;
   
### Running the Pipeline

```bash
chmod +x runMaxEnt.sh aux/*.sh aux/*.py
nohup ./runMaxEnt.sh & 
```

&nbsp;

The script will:
1. Run MaxEnt for M=0 through M=M_max (set in the script)
2. Collect results into summary files
3. Average decay width functions over different M ranges
4. Display final results with uncertainties

&nbsp;

## Output Files

&nbsp;

### Summary Data Files

| File | Description |
|------|-------------|
| `Gm.dat` | Ground state Gm(0) values for each M |
| `GmEres.dat` | Decay width Gamma(Er) at resonance energy |
| `DlEres.dat` | Level shift Delta(Er) at resonance energy |
| `EdEr.dat` | Discrete state (Ed) and resonance (Er) energies |
| `Dmu.dat` | Moment differences (input vs reconstructed from MaxEnt) |
| `ent.dat` | Final values for entropies S_x and S_E |
| `etas.2d` | Lagrange multipliers (eta values) from all runs; use the `eta.plt` script from the `gnuplot` directory to visualize the data |

&nbsp;

### Averaged Decay width Functions

| File | M Range | Use Case |
|------|---------|----------|
| `func.aver.03_20` | M=3-20 | Maximum statistics |
| `func.aver.05_15` | M=5-15 | Recommended range |
| `func.aver.05_20` | M=5-20 | High-moment range |

Each averaged file contains:
- Column 1: Energy E
- Column 2: Average Gamma(E)
- Column 3: Standard deviation
- Column 4: Standard error
- Footer: Gamma(0) ± uncertainty [meV]
- Footer: Lifetime tau ± uncertainty [fs]

&nbsp;

### Run Directories

Each `run##` directory contains:
- `ent.out` - Complete MaxEnt output
- `fort.3000` - Decay width function Gamma(E) in meV vs energy in eV
- `fort.9006` - spectral moments convergence data

&nbsp;

## Script Details

&nbsp;

### collect.sh

**Extracts quantities from output files**

Searches for specific patterns in `ent.out`:
- `Gm(0)_ent` - Decay width at the discrete state energy
- `Gamma(Er)` - Decay Width at the resonance energy
- `Delta(Er)` - Level shift at the resonance energy
- `dmu2` - convergence of the spectral moments
- `entropy_ent` - Entropies in X and E variable

&nbsp;

### inp_subs.py

**Modifies input file parameters**

Usage: `inp_subs.py <file> <variable> <value>`

Safely replaces variable assignments in Fortran namelist style:
```
NT = 00   →   NT = 05
```

Creates `<file>.new` with modifications.

&nbsp;

### get_max_dmu.py

**Analyzes moment difference convergence**

Processes `fort.9006` files containing input and reconstructed spectral moments
- Finds the moment with maximum error `dmu_max` and the corresponding `k`
- Computes average error across all k's
- Outputs: `<k> <dmu_max> <dmu_avg>`

&nbsp;

### maxentAverage.py

**Computes ensemble-averaged decay width functions**

Usage: `maxentAverage.py func.all <M_low> <M_high>`

Methodology:
1. Extracts Gamma columns for M ∈ [M_low, M_high]
2. Computes mean and standard deviation at each energy
3. Interpolates to find Gamma(0)
4. Converts to lifetime: τ = ℏ/Γ = 658.21 fs·meV / Γ[meV]
