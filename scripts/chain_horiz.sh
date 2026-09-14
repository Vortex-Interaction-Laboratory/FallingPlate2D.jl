#!/bin/bash
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:$PATH"
for D in 1e-5 1e-4; do
  [ -f "ckpt/flap_$D.jld2" ] || { echo "MISSING ckpt/flap_$D.jld2"; exit 1; }
  for A in 0.5 1.5 2.5; do
    for o in pos neg; do
      nm="h_${D}_A${A}_${o}"
      if [ -s "out/$nm.csv" ] && [ "$(tail -1 out/$nm.csv | cut -d, -f1 | cut -d. -f1)" = "40" ]; then
        echo "skip $nm (complete)"; continue
      fi
      echo "=== $nm start $(date +%T)"
      julia -t 4 --project=. scripts/matrix_horiz.jl "$D" "$A" "$o" > "logs/$nm.log" 2>&1
      echo "=== $nm end   $(date +%T) exit=$?"
    done
  done
done
echo "HORIZ MATRIX DONE $(date +%T)"
