#!/bin/bash
cargo build --release --example cachemisses

PROJECT_DIR=$(dirname "$(cargo locate-project --workspace --message-format plain)")
TARGET_DIR="$PROJECT_DIR/target/release/examples"
cd $TARGET_DIR

echo "n reps presort f32 Ir Dr Dw I1mr D1mr D1mw ILmr DLmr DLmw"

exponents="2 3 4 5 6 7 8"


# $1: [true|false] ... pre-sort random test data
# $2: [true|false] ... single-precision data & computation (i.e. [f32|f64])
presort=${1:-false}
f32=${2:-false}


for i in $exponents
do
  n=$((10**i))
  # repeats each CellGrid::new() $reps times to get comparable runtimes
  # (each command then takes ~30-60s)
  # reps=$((10**(7-i)))
  reps=1

  values=$(valgrind --tool=callgrind --callgrind-out-file=./cachemisses_$n.out \
    --cache-sim=yes --instr-atstart=no ./cachemisses \
    $n $reps $presort $f32 \
    2>&1 >/dev/null \
    | grep "Collected" \
    | grep -oE '( [[:digit:]]+){9}')
  echo "$n $reps $presort $f32$values"
done
