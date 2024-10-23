#!/bin/bash

# Usage: ./generate_topologies.sh input/ligA.mol2 input/ligB.mol2 output/data

LIG_A_PATH=$1
LIG_B_PATH=$2
OUTPUT_DIR=$3

if [ -z "$OUTPUT_DIR" ]; then
  OUTPUT_DIR="output/data"
fi

mkdir -p $OUTPUT_DIR

process_ligand() {
    LIGAND_PATH=$1
    BASE_NAME=$2

    acpype -i "$LIGAND_PATH" -b "$BASE_NAME" -c bcc -a gaff2 -o gmx -w
    mv "${BASE_NAME}.acpype/${BASE_NAME}_GMX.itp" "${OUTPUT_DIR}/${BASE_NAME}.itp"
    rm -r "${BASE_NAME}.acpype"
}

process_ligand "$LIG_A_PATH" "ligA"
process_ligand "$LIG_B_PATH" "ligB"

echo "All ligand topologies have been successfully generated."
