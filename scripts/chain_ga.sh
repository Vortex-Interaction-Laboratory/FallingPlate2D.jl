#!/bin/bash
cd "$(dirname "$0")/.."; export PATH="/opt/homebrew/bin:$PATH"
julia -t 2 --project=. scripts/matrix_ga.jl lap 50  > logs/galap_Ga50.log  2>&1 &
A=$!; julia -t 2 --project=. scripts/matrix_ga.jl lap 100 > logs/galap_Ga100.log 2>&1 &
B=$!; wait $A; RA=$?; wait $B; RB=$?
echo "GA LAPS DONE $(date +%T) exit=$RA,$RB"
[ $RA -eq 0 ] && [ $RB -eq 0 ] || exit 1
run(){ nm="g_Ga$1_$2${3:+_$3}"; [ "$2" = const ] && nm="g_Ga$1_const"
  [ -s "out/$nm.csv" ] && [ "$(tail -1 out/$nm.csv|cut -d, -f1|cut -d. -f1)" = 40 ] && { echo "skip $nm"; return; }
  echo "=== $nm start $(date +%T)"
  julia -t 4 --project=. scripts/matrix_ga.jl run "$1" "$2" "$3" > "logs/$nm.log" 2>&1
  echo "=== $nm end   $(date +%T) exit=$?"; }
for G in 50 100; do
  run $G const
  run $G v down; run $G v up
  run $G h pos;  run $G h neg
done
echo "GA MATRIX DONE $(date +%T)"
