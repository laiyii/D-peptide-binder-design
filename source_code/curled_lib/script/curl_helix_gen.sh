#!/bin/bash

output_name=""
outdir=""
len=""
phi0=""
delphi=""
psi0=""
delpsi=""
phase=""

while [[ "$1" != "--" ]]; do
  case $1 in
    -o)
      output_name=$2
      shift 2
      ;;
    -outdir)
      outdir=$2
      shift 2
      ;;
    -len)
      len=$2
      shift 2
      ;;
    -phi0)
      phi0=$2
      shift 2
      ;;
    -delphi)
      delphi=$2
      shift 2
      ;;
    -psi0)
      psi0=$2
      shift 2
      ;;
    -delpsi)
      delpsi=$2
      shift 2
      ;;
    -phase)
      phase=$2
      shift 2
      ;;
    *)
      echo "invalid options: $1" >&2
      exit 1
      ;;
  esac
done

shift

if [ -z "$output_name" ] || [ -z "$outdir" ] || [ -z "$len" ] || [ -z "$phi0" ] || [ -z "$delphi" ] || [ -z "$psi0" ] || [ -z "$delpsi" ] || [ -z "$phase" ]; then
  echo "All options -o, -outdir, -len, -phi0, -delphi, -psi0, -delpsi, and -phase must be specified" >&2
  exit 1
fi

mkdir -p "$outdir"

"$HSD/curled_lib/script/PhiPsi2Helix" "$outdir/$output_name" "$len" "$phi0" "$delphi" "$psi0" "$delpsi" "$phase"