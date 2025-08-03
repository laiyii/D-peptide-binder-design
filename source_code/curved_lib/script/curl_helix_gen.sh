#!/bin/bash

#===============================================================================
# CURL HELIX GENERATOR - BATCH PROCESSING SCRIPT
#===============================================================================
#
# DESCRIPTION:
#   This script generates multiple helical polyALA structures by reading 
#   dihedral angle parameters from a CSV file and creating all possible 
#   parameter combinations. It calls the PhiPsi2Helix program to generate 
#   PDB files for each unique combination of phi, psi, and phase angles.
#
# USAGE:
#   ./curl_helix_gen.sh -outdir <output_directory> -len <length> -paramlist <csv_file>
#
# PARAMETERS:
#   -outdir     : Output directory where PDB files will be saved
#   -len        : Length of the helix (number of residues)  
#   -paramlist  : CSV file containing parameter combinations
#
# CSV FILE FORMAT:
#   The CSV file must have the following header and structure:
#   
#   phi0,delphi,psi0,delpsi,phase
#   -62.0,-3.0,-39.0,-1.0,-60.0
#   -60.0,-2.0,-40.0,-1.5,-50.0
#   -65.0,-4.0,-35.0,-0.5,-70.0
#
#   Where:
#   - phi0    : Initial phi dihedral angle (degrees)
#   - delphi  : Phi angle increment per residue (degrees)
#   - psi0    : Initial psi dihedral angle (degrees) 
#   - delpsi  : Psi angle increment per residue (degrees)
#   - phase   : Phase angle for helical twist (degrees)
#
# OUTPUT:
#   - PDB files named: H_<len>_<phi0>_<delphi>_<psi0>_<delpsi>_<phase>.pdb
#   - Decimal points are replaced with underscores
#   - Negative signs are replaced with 'm'
#   - Example: H_21_m62_0_m3_0_m39_0_m1_0_m60_0.pdb
#            (from phi0=-62.0, delphi=-3.0, psi0=-39.0, delpsi=-1.0, phase=-60.0)
#
# NOTE THAT:
#   - Only when BOTH delphi=0 AND delpsi=0, the phase parameter becomes meaningless
#   - In this case, the script automatically uses only phase=0
#   - If only one of them is 0, all phase values will still be traversed
#   - This ensures proper sampling while avoiding truly redundant calculations
#
#
# REQUIREMENTS:
#   - PhiPsi2Helix executable must be available at $DPEP/curled_lib/script/
#   - CSV file with proper format and valid numerical parameters
#   - Write permissions for the output directory
#
# AUTHOR: Laiyi Feng, Changsheng Zhang
# VERSION: 1.0
# DATE: 2021-12-10
#
#===============================================================================

outdir=""
len=""
paramlist=""

while [[ $# -gt 0 && "$1" != "--" ]]; do
  case $1 in
    -outdir)
      if [[ -n "$2" && "$2" != -* ]]; then
        outdir=$2
        shift 2
      else
        echo "Error: -outdir requires a directory path" >&2
        exit 1
      fi
      ;;
    -len)
      if [[ -n "$2" && "$2" != -* ]]; then
        len=$2
        shift 2
      else
        echo "Error: -len requires a number" >&2
        exit 1
      fi
      ;;
    -paramlist)
      if [[ -n "$2" && "$2" != -* ]]; then
        paramlist=$2
        shift 2
      else
        echo "Error: -paramlist requires a CSV file path" >&2
        exit 1
      fi
      ;;
    -h|--help)
      echo "Usage: $0 -outdir <output_directory> -len <length> -paramlist <csv_file>"
      echo ""
      echo "Options:"
      echo "  -outdir     Output directory where PDB files will be saved"
      echo "  -len        Length of the helix (number of residues)"
      echo "  -paramlist  CSV file containing parameter combinations"
      echo "  -h, --help  Show this help message"
      exit 0
      ;;
    *)
      echo "Error: Unknown option '$1'" >&2
      echo "Usage: $0 -outdir <output_directory> -len <length> -paramlist <csv_file>"
      echo "Use -h or --help for more information"
      exit 1
      ;;
  esac
done

# Check required parameters
if [ -z "$outdir" ] || [ -z "$len" ] || [ -z "$paramlist" ]; then
  echo "All options -outdir, -len, and -paramlist must be specified" >&2
  echo "Usage: $0 -outdir <output_directory> -len <length> -paramlist <csv_file>"
  exit 1
fi

# Check if CSV file exists
if [ ! -f "$paramlist" ]; then
  echo "Parameter file '$paramlist' not found" >&2
  exit 1
fi

mkdir -p "$outdir"

# Read CSV file and extract unique values for each parameter
echo "Reading parameters from $paramlist..."

# Check if the CSV has a header and extract column indices
header=$(head -n 1 "$paramlist")
echo "CSV header: $header"

# Validate CSV format and data
echo "Validating CSV data..."
line_count=$(wc -l < "$paramlist")
echo "Total lines in CSV: $line_count (including header)"

# Check for common CSV issues
echo "Sample data lines (first 3 data rows):"
tail -n +2 "$paramlist" | head -3 | while IFS= read -r line; do
  echo "  '$line'"
  # Check for inconsistent field counts
  field_count=$(echo "$line" | tr ',' '\n' | wc -l)
  if [ "$field_count" -ne 5 ]; then
    echo "    WARNING: Expected 5 fields, found $field_count"
  fi
done

# Extract unique values for each parameter column
phi0_values=($(tail -n +2 "$paramlist" | cut -d',' -f1 | sed 's/[[:space:]]//g' | sort -n | uniq))
delphi_values=($(tail -n +2 "$paramlist" | cut -d',' -f2 | sed 's/[[:space:]]//g' | sort -n | uniq))
psi0_values=($(tail -n +2 "$paramlist" | cut -d',' -f3 | sed 's/[[:space:]]//g' | sort -n | uniq))
delpsi_values=($(tail -n +2 "$paramlist" | cut -d',' -f4 | sed 's/[[:space:]]//g' | sort -n | uniq))
phase_values=($(tail -n +2 "$paramlist" | cut -d',' -f5 | sed 's/[[:space:]]//g' | sort -n | uniq))

# Display parameter ranges
echo "Parameter ranges found:"
echo "  phi0: ${phi0_values[@]} (${#phi0_values[@]} values)"
echo "  delphi: ${delphi_values[@]} (${#delphi_values[@]} values)"
echo "  psi0: ${psi0_values[@]} (${#psi0_values[@]} values)"
echo "  delpsi: ${delpsi_values[@]} (${#delpsi_values[@]} values)"
echo "  phase: ${phase_values[@]} (${#phase_values[@]} values)"

total_combinations=$((${#phi0_values[@]} * ${#delphi_values[@]} * ${#psi0_values[@]} * ${#delpsi_values[@]} * ${#phase_values[@]}))
echo "Total combinations to generate: $total_combinations (before optimization)"

# Calculate optimized combinations
optimized_combinations=0
for phi0 in "${phi0_values[@]}"; do
  for delphi in "${delphi_values[@]}"; do
    for psi0 in "${psi0_values[@]}"; do
      for delpsi in "${delpsi_values[@]}"; do
        # Check if we need to use only phase=0 (only when BOTH are zero)
        if [[ "$delphi" == "0" || "$delphi" == "0.0" ]] && [[ "$delpsi" == "0" || "$delpsi" == "0.0" ]]; then
          optimized_combinations=$((optimized_combinations + 1))
        else
          optimized_combinations=$((optimized_combinations + ${#phase_values[@]}))
        fi
      done
    done
  done
done

echo "Optimized combinations to generate: $optimized_combinations"
total_combinations=$optimized_combinations

# Generate all combinations
counter=0
failed_counter=0

for phi0 in "${phi0_values[@]}"; do
  for delphi in "${delphi_values[@]}"; do
    for psi0 in "${psi0_values[@]}"; do
      for delpsi in "${delpsi_values[@]}"; do
        
        # Optimization logic for phase parameter:
        # - If both delphi=0 AND delpsi=0: only use phase=0
        # - If only one of them is 0: use all phase values (full traversal)
        # - If both are non-zero: use all phase values (normal case)
        if [[ "$delphi" == "0" || "$delphi" == "0.0" ]] && [[ "$delpsi" == "0" || "$delpsi" == "0.0" ]]; then
          # Both increments are 0, use only phase=0
          phase_list=(0)
          optimization_note="Both delphi and delpsi are zero"
        else
          # Either one or both increments are non-zero, use all phase values
          phase_list=("${phase_values[@]}")
          optimization_note=""
        fi
        
        for phase in "${phase_list[@]}"; do
          # Generate output filename
          output_name="H_${len}_${phi0}_${delphi}_${psi0}_${delpsi}_${phase}"
          
          # Clean and format the filename:
          # 1. Remove all whitespace characters
          # 2. Replace negative signs with 'm'
          # 3. Replace decimal points with underscores
          # 4. Add .pdb extension
          output_name=$(echo "$output_name" | sed 's/[[:space:]]//g' | sed 's/-/m/g' | sed 's/\./_/g')
          output_name="${output_name}.pdb"
          
          # Debug: Show raw parameters to help identify issues
          counter=$((counter + 1))
          echo "[$counter/$total_combinations] Generating: $output_name"
          echo "  Raw parameters: phi0='$phi0', delphi='$delphi', psi0='$psi0', delpsi='$delpsi', phase='$phase'"
          
          # Show optimization info
          if [[ -n "$optimization_note" ]]; then
            echo "  Note: Using phase=0 only ($optimization_note)"
          fi
          
          # Call the PhiPsi2Helix program
          if "$DPEP/curled_lib/script/PhiPsi2Helix" "$outdir/$output_name" "$len" "$phi0" "$delphi" "$psi0" "$delpsi" "$phase"; then
            echo "  ✓ Success: $output_name"
          else
            echo "  ✗ Failed: $output_name"
            failed_counter=$((failed_counter + 1))
          fi
          echo ""
        done
      done
    done
  done
done

echo "========================================="
echo "Batch generation completed!"
echo "Total structures requested: $total_combinations"
echo "Successfully generated: $((counter - failed_counter))"
echo "Failed: $failed_counter"
echo "Output directory: $outdir"
echo "========================================="