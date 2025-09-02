#!/bin/bash
RUSTFLAGS="-C target-cpu=native" cargo build --package psssh --release --example psssh

PROJECT_DIR=$(dirname "$(cargo locate-project --workspace --message-format plain)")
TARGET_DIR="$PROJECT_DIR/target/release/examples"

echo "PDB, n, vaabb, cutoff, queries, ns"

AXIS_LENGTH=64
cutoffs="1 2 5 10"

for c in $cutoffs
do
  for pdb in "$@"
  do
    output=$("$TARGET_DIR/psssh" eval -c $c -l $AXIS_LENGTH $pdb)
    echo $output
  done
done
