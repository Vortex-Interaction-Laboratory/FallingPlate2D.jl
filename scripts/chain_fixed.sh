#!/bin/bash
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:$PATH"
julia -t 2 --project=. scripts/matrix_fixed.jl lap 1e-5 > logs/flap_1e-5.log 2>&1 &
A=$!
julia -t 2 --project=. scripts/matrix_fixed.jl lap 1e-4 > logs/flap_1e-4.log 2>&1 &
B=$!
wait $A; RA=$?; wait $B; RB=$?
echo "FIXED LAPS DONE $(date +%T) exit=$RA,$RB"
[ $RA -eq 0 ] && [ $RB -eq 0 ] || { echo "lap failed"; exit 1; }
run() {
  local nm; if [ "$2" = "const" ]; then nm="f_$1_const"; else nm="f_$1_A$2_$3"; fi
  [ -s "out/$nm.csv" ] && { echo "skip $nm"; return; }
  echo "=== $nm start $(date +%T)"
  julia -t 4 --project=. scripts/matrix_fixed.jl run "$1" "$2" "$3" > "logs/$nm.log" 2>&1
  echo "=== $nm end   $(date +%T) exit=$?"
}
for D in 1e-5 1e-4; do
  run "$D" const
  for A in 0.5 1.5 2.5; do run "$D" "$A" down; run "$D" "$A" up; done
done
echo "FIXED MATRIX DONE $(date +%T)"
