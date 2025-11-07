#!/bin/bash
cargo build --release --example lmp-data

# $1: [true|false] ... use LAMMPS (default) or CellListMap.jl
use_lammps=${1:-true}

TMP_DIR="/tmp/zelll_benches"
>&2 echo "creating $TMP_DIR for input data"
mkdir $TMP_DIR

LAMMPS=lmp
PROJECT_DIR=$(dirname "$(cargo locate-project --workspace --message-format plain)")
LMP_COMMANDS="$PROJECT_DIR/more_benches/in.zelllbench.txt"
cd $PROJECT_DIR

# 10⁸ requires large amounts of RAM (>20GB)
exponents=(2 3 4 5 6 7 8)
repeats=(100000 10000 1000 100 10 1 1)
# use this instead if RAM is limited
#exponents=(2 3 4 5 6)
#repeats=(100000 10000 1000 100 10)

echo "n reps energy runtime memory tool"

for i in "${!exponents[@]}"
do
  exponent=${exponents[$i]}
  repeat=${repeats[$i]}
  n=$((10**exponent))

  TMP_FILE="$TMP_DIR/${n}atomsinabox.txt"
  >&2 echo "creating data: $TMP_FILE"
  cargo run --release --example lmp-data -- $n > $TMP_FILE

  if ${use_lammps}; then
    dump=$($LAMMPS -in $LMP_COMMANDS -var data $TMP_FILE -var repeat $repeat)

    memory=$(echo $dump \
      | grep -oE '| ([[:digit:]]+.[[:digit:]]+)(e[+,-]?[[:digit:]]+)? Mbytes' \
      | grep -oE '([[:digit:]]+.[[:digit:]]+)(e[+,-]?[[:digit:]]+)?')

    time=$(echo $dump \
      | grep -oE 'Loop time of ([[:digit:]]+.[[:digit:]]+) on' \
      | grep -oE '([[:digit:]]+.[[:digit:]]+)')

    e_pot=$(echo $dump \
      | grep -oE '( - \[0, 0, 0, [[:digit:]]+.[[:digit:]]+)(e[+,-]?[[:digit:]]+)?,' \
      | grep -oE '([[:digit:]]+.[[:digit:]]+)(e[+,-]?[[:digit:]]+)?')

    echo $n $repeat $e_pot $time $memory "\"LAMMPS\""
  else
    echo $(julia --project=$PROJECT_DIR/more_benches/ $PROJECT_DIR/more_benches/celllistmap.jl $TMP_FILE)
  fi
done

>&2 echo "deleting $TMP_DIR and its content"
rm -r $TMP_DIR