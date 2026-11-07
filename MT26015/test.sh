#!/bin/bash
# Checks every CSV is identical to the sequential one. Usage: ./test.sh
set -e
g++ -O2 -std=c++17 -pthread -o seq seq.cpp && g++ -O2 -std=c++17 -pthread -o par par.cpp
A="--N 150 --days 120 --seed 7 --beta 0.12 --v 0.1 --m 0.9 --I0 20 --fixed 1"
./seq $A --out ref.csv 2>/dev/null
fail=0; n=0
for T in 1 2 3 4 7 8 16; do for s in strips blocks interleaved; do
  for cfg in "thread $T" "pool $T" "pool $((3*T+1))"; do
    set -- $cfg; ./par $A --mode $1 --threads $T --chunks $2 --split $s --out t.csv 2>/dev/null
    n=$((n+1)); cmp -s ref.csv t.csv || { echo "MISMATCH: mode=$1 T=$T C=$2 split=$s"; fail=1; }
  done; done; done
./par $A --mode seq --out t.csv 2>/dev/null; cmp -s ref.csv t.csv || { echo "MISMATCH par-seq"; fail=1; }
[ $fail = 0 ] && echo "ALL $((n+1)) runs identical to sequential" || exit 1