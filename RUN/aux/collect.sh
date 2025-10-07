#!/bin/bash
#
# collect.sh - Extract and organize MaxEnt results from run directories
#
# Purpose: Collects key quantities (Gm, Gamma, Delta, entropy, etc.) from 
#          MaxEnt output files and organizes them into summary data files
#
# Output files:
#   - Gm.dat: Gm(0)_ent values
#   - GmEres.dat: Gamma(Er) values at resonance energy
#   - DlEres.dat: Delta(Er) values at resonance energy
#   - EdEr.dat: Ed and Er energies
#   - Dmu.dat: Moment differences (input vs reconstructed moments)
#   - ent.dat: Entropy values
#

set -e  # Exit on error
set -u  # Exit on undefined variable

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# Extract Gm(0)_ent - Ground state values
# ============================================================================

echo "Extracting Gm(0)_ent values..."
{
    for rundir in run??; do
        if [[ ! -d "$rundir" ]]; then
            continue
        fi
        
        # Extract run number (remove "run" prefix)
        run_num="${rundir#run}"
        
        # Extract last occurrence of Gm(0)_ent from output
        if [[ -f "${rundir}/ent.out" ]]; then
            gm=$(grep "Gm(0)_ent" "${rundir}/ent.out" | tail -n 1 || echo "")
            if [[ -n "$gm" ]]; then
                echo "$run_num $gm"
            else
                echo "Warning: No Gm(0)_ent found in ${rundir}/ent.out" >&2
            fi
        else
            echo "Warning: ${rundir}/ent.out not found" >&2
        fi
    done
} > Gm.dat

# ============================================================================
# Extract Gamma(Er) - Width at resonance energy
# ============================================================================

echo "Extracting Gamma(Er) values..."
{
    for rundir in run??; do
        if [[ ! -d "$rundir" ]]; then
            continue
        fi
        
        run_num="${rundir#run}"
        
        if [[ -f "${rundir}/ent.out" ]]; then
            gm=$(grep "Gamma(Er)" "${rundir}/ent.out" | tail -n 1 || echo "")
            if [[ -n "$gm" ]]; then
                echo "$run_num $gm"
            else
                echo "Warning: No Gamma(Er) found in ${rundir}/ent.out" >&2
            fi
        fi
    done
} > GmEres.dat

# ============================================================================
# Extract Delta(Er) - Level shift at resonance energy
# ============================================================================

echo "Extracting Delta(Er) values..."
{
    for rundir in run??; do
        if [[ ! -d "$rundir" ]]; then
            continue
        fi
        
        run_num="${rundir#run}"
        
        if [[ -f "${rundir}/ent.out" ]]; then
            gm=$(grep "Delta(Er)" "${rundir}/ent.out" | tail -n 1 || echo "")
            if [[ -n "$gm" ]]; then
                echo "$run_num $gm"
            else
                echo "Warning: No Delta(Er) found in ${rundir}/ent.out" >&2
            fi
        fi
    done
} > DlEres.dat

# ============================================================================
# Extract Moment Differences (dmu = difference between input and recovered)
# ============================================================================

echo "Extracting moment differences (input vs reconstructed)..."

# Remove old file if exists
rm -f EdEr.dat

# Create header for Dmu.dat
# Here, M is the highest moment in the run, D_init/D_final are the initial and
# final values of the minimized D(lambda) function, k_w and dmu_w is the index
# and error of the worst reconstructed spectral moment, and dmu_aver is the
# average error across all spectral momenst
cat > Dmu.dat << 'EOF'
# M  D_final              D_init               k_w  dmu_w     dmu_aver
==========================================================================
EOF

for rundir in run??; do
    if [[ ! -d "$rundir" ]]; then
        continue
    fi
    
    run_num="${rundir#run}"
    
    if [[ -f "${rundir}/ent.out" ]]; then
        # Extract final dmu2 value (moment difference at convergence)
        dmu=$(grep "dmu2" "${rundir}/ent.out" | tail -n 1 | awk '{print $NF}' || echo "N/A")
        
        # Extract initial dmu2 value (moment difference at start)
        dmu0=$(grep "dmu2" "${rundir}/ent.out" | head -n 1 | awk '{print $NF}' || echo "N/A")
        
        # Extract maximum moment difference statistics using Python script
        if [[ -f "${rundir}/fort.9006" && -x "${SCRIPT_DIR}/get_max_dmu.py" ]]; then
            dmumax=$("${SCRIPT_DIR}/get_max_dmu.py" "${rundir}/fort.9006" || echo "N/A N/A N/A")
        else
            dmumax="N/A N/A N/A"
        fi
        
        echo "$run_num $dmu $dmu0 $dmumax"
        
        # Extract Ed and Er energies
        Ed=$(grep "Ed \[eV\]" "${rundir}/ent.out" | tail -n 1 | awk '{print $NF}' || echo "N/A")
        Er=$(grep "Er \[eV\]" "${rundir}/ent.out" | tail -n 1 | awk '{print $NF}' || echo "N/A")
        
        echo "$run_num $Ed $Er" >> EdEr.dat
    else
        echo "Warning: ${rundir}/ent.out not found" >&2
    fi
done >> Dmu.dat

# ============================================================================
# Extract Entropy Values
# ============================================================================

echo "Extracting entropy values..."

# Create header for ent.dat
cat > ent.dat << 'EOF'
M                S_x[Gamma]        S_E[Gamma]
=============================================
EOF

for rundir in run??; do
    if [[ ! -d "$rundir" ]]; then
        continue
    fi
    
    run_num="${rundir#run}"
    
    if [[ -f "${rundir}/ent.out" ]]; then
        ent=$(grep "entropy_ent" "${rundir}/ent.out" | tail -n 1 || echo "")
        if [[ -n "$ent" ]]; then
            echo "$run_num $ent"
        else
            echo "Warning: No entropy_ent found in ${rundir}/ent.out" >&2
        fi
    fi
done >> ent.dat

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "Collection complete. Generated files:"
echo "  - Gm.dat: Ground state Gm(0) values"
echo "  - GmEres.dat: Gamma at resonance energy"
echo "  - DlEres.dat: Delta at resonance energy"
echo "  - EdEr.dat: Ed and Er energies"
echo "  - Dmu.dat: Moment differences (input vs reconstructed)"
echo "  - ent.dat: Entropy values"
echo ""
