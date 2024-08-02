#!/bin/bash

input_file=""
output_file=""

while getopts "i:o:" opt; do
  case $opt in
    i)
      input_file=$OPTARG
      ;;
    o)
      output_file=$OPTARG
      ;;
    \?)
      echo "invalid options: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "-$OPTARG needs an input" >&2
      exit 1
      ;;
  esac
done

if [ -z "$input_file" ]; then
  echo "must use -i to specify the input file" >&2
  exit 1
fi

filename=$(basename "$input_file" .pdb)

if [ -z "$output_file" ]; then
  output_file="${filename}_mirror.pdb"
fi


echo "input file: $input_file"
echo "output file: $output_file"

python3 $HSD/docking/mirror_target/mirror_protein.py -i "$input_file" -o "$output_file" && echo "flipped pdb file generated"