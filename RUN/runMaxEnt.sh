#!/bin/bash
#
# runMaxEnt.sh - Driver script for Maximum Entropy analysis
#
# Purpose: Runs MaxEnt calculations for M=0-maxM, collecting and averaging results
# Dependencies: entropy.exe, readEta.exe, auxiliary Python scripts
# 

maxM=15
inpfile=$(grep -i INPFILE entropy.input | awk -F'"' '{ print $2 }')
if [ -z ${inpfile} ]
then
  echo "Error: cannot extract input data file name from entropy.input"
  exit 1
fi


# ============================================================================
# Setup and Environment
# ============================================================================

# Set Intel compiler environment, OMP_NUM_THREADS, stack limits, etc.
# Note: Source envir.sh before setting -u because Intel scripts may reference $1
if [[ -f ./envir.sh ]]; then
    source ./envir.sh
else
    echo "Warning: ./envir.sh not found - environment may not be properly configured"
fi

# Get absolute path to parent directory (where executables are located)
# This remains fixed even when we change directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exedir="$(cd "${SCRIPT_DIR}/.." && pwd)"
executable="${exedir}/entropy.exe"
auxdir="${SCRIPT_DIR}/aux"

# Verify executables exist
if [[ ! -x "$executable" ]]; then
    echo "Error: Executable not found or not executable: $executable"
    exit 1
fi

if [[ ! -x "${exedir}/readEta.exe" ]]; then
    echo "Error: readEta.exe not found or not executable: ${exedir}/readEta.exe"
    exit 1
fi

# Save original working directory
WORK_DIR="$PWD"

# ============================================================================
# Run MaxEnt for M=0-maxM
# ============================================================================

echo "Starting MaxEnt calculations for M=0 to ${maxM}..."
echo "Working directory: ${WORK_DIR}"
echo "Executable: ${executable}"
echo ""
echo "Configuration: "
cat entropy.input
echo ""

for nt in $(seq -w 00 ${maxM}); do
    rundir="run${nt}"
    echo "Processing M=${nt} in ${rundir}..."
    
    # Create run directory and copy input files
    mkdir -p "${rundir}"
    
    # Copy required input files
    if [[ ! -f entropy.input ]]; then
        echo "Error: entropy.input not found in ${WORK_DIR}"
        exit 1
    fi
    cp entropy.input "${rundir}/"
    
    if [[ ! -f "${inpfile}" ]]; then
        echo "Error: ${inpfile} not found in ${WORK_DIR}"
        exit 1
    fi
    cp ${inpfile} "${rundir}/"
    
    # Copy eta.last if it exists (optional - code runs with or without it)
    # Typically it contains Lagrange multiplicators found in a previous run with
    # M=NT-1
    if [[ -f eta.last ]]; then
        cp eta.last "${rundir}/"
    fi
    
    cd "${rundir}"
    
    # Modify input file to set the actual NT parameter (M)
    if [[ -x "${auxdir}/inp_subs.py" ]]; then
        "${auxdir}/inp_subs.py" entropy.input NT "${nt}"
        mv -f entropy.input.new entropy.input
    else
        echo "Error: inp_subs.py not found or not executable"
        cd "${WORK_DIR}"
        exit 1
    fi
    
    # Run MaxEnt calculation with timing
    echo "  Running entropy.exe..."
    if /usr/bin/time -v "$executable" < entropy.input > ent.out 2>&1; then
        echo "  Completed successfully"
        
        # Copy eta.last back to parent for next run (overwrite without prompt)
        if [[ -f eta.last ]]; then
            cp -f eta.last ../.
        fi
    else
        echo "  ERROR: entropy.exe failed for M=${nt}"
        echo "  Check ${rundir}/ent.out for details"
        cd "${WORK_DIR}"
        exit 1
    fi
    
    cd "${WORK_DIR}"
done

echo ""
echo "All MaxEnt calculations completed."
echo ""

# ============================================================================
# Collect and Display Results
# ============================================================================

echo "Collecting results from all runs..."

if [[ -x "${auxdir}/collect.sh" ]]; then
    "${auxdir}/collect.sh"
else
    echo "Error: collect.sh not found or not executable"
    exit 1
fi

# Display summary of collected data files
echo ""
echo "Summary of results:"
for datfile in *.dat; do
    if [[ -f "$datfile" ]]; then
        echo "--- ${datfile} (last 5 lines) ---"
        tail -n 5 "${datfile}"
        echo ""
    fi
done

# ============================================================================
# Extract Lagrange Multipliers (Lambda Values)
# ============================================================================

echo "Extracting Lagrange multipliers from individual runs..."

rm -f etas.2d
for rundir in run??; do
    if [[ -d "$rundir" ]]; then
        cd "${rundir}"
        "${exedir}/readEta.exe" 35 >> "${WORK_DIR}/etas.2d"
        cd "${WORK_DIR}"
    fi
done

echo "Lagrange multipliers saved to etas.2d"
echo ""

# ============================================================================
# Average Gamma(E) Over Different M Ranges
# ============================================================================

echo "Computing averaged Gamma(E) for various M ranges..."

MEaverpy="${auxdir}/maxentAverage.py"

if [[ ! -x "$MEaverpy" ]]; then
    echo "Error: maxentAverage.py not found or not executable"
    exit 1
fi

# Collect all fort.3000 files into one
if paste run*/fort.3000 >| func.all 2>/dev/null; then
    echo "Created func.all from fort.3000 files"
else
    echo "Warning: Could not create func.all - some runs may be missing fort.3000"
fi

# Define M ranges for averaging
declare -a ranges=(
    "03 15"
    "04 15"
    "05 15"
    "03 20"
    "05 20"
    "04 20"
)

echo ""
for range in "${ranges[@]}"; do
    read -r low high <<< "$range"
    outfile="func.aver.${low}_${high}"
    echo "Averaging M=${low} to M=${high} -> ${outfile}"
    "$MEaverpy" func.all "$low" "$high" > "$outfile"
done

# ============================================================================
# Display Final Averaged Results
# ============================================================================

echo ""
echo "=========================================="
echo "Final averaged Gamma_loc with std dev:"
echo "=========================================="
tail -n 2 func.aver.*

echo ""
echo "MaxEnt analysis pipeline completed successfully!"
