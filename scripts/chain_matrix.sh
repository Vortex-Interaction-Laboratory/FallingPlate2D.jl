#!/bin/bash
# Wait for the Da=1e-5 and Da=1e-3 laps, then run the full matrix sequentially.
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:$PATH"

for i in $(seq 1 720); do
  grep -q "LAPS DONE" logs/chain_laps.log 2>/dev/null && break
  sleep 10
done
grep -q "LAPS DONE" logs/chain_laps.log || { echo "laps never finished"; exit 1; }
for D in 1e-5 1e-3; do
  [ -f "ckpt/mlap_$D.jld2" ] || { echo "MISSING ckpt/mlap_$D.jld2"; exit 1; }
done
echo "laps ready $(date +%T); starting matrix"

run() {  # $1=Da  $2=A|const  $3=order
  local nm
  if [ "$2" = "const" ]; then nm="m_$1_const"; else nm="m_$1_A$2_$3"; fi
  if [ -s "out/$nm.csv" ]; then echo "skip $nm (exists)"; return; fi
  echo "=== $nm  start $(date +%T)"
  julia -t 4 --project=. scripts/matrix.jl run "$1" "$2" "$3" > "logs/$nm.log" 2>&1
  echo "=== $nm  end   $(date +%T)  exit=$?"
}

for D in 1e-5 1e-4 1e-3; do
  run "$D" const
  for A in 0.5 1.5 2.5; do
    run "$D" "$A" down
    run "$D" "$A" up
  done
done
echo "MATRIX DONE $(date +%T)"
