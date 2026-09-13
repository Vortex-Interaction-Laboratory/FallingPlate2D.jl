#!/bin/bash
# Wait for lap1, then branch end1 and end2 in parallel, then plot.
cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:$PATH"

echo "waiting for lap1 ..."
for i in $(seq 1 360); do
  grep -q "saved ckpt/lap1.jld2" logs/lap1.log 2>/dev/null && break
  grep -q "aborted:" logs/lap1.log 2>/dev/null && { echo "lap1 ABORTED"; exit 1; }
  sleep 10
done
[ -f ckpt/lap1.jld2 ] || { echo "no checkpoint; giving up"; exit 1; }
echo "lap1 done at $(date +%T)"

# two independent branches from the same state, in parallel, 2 threads each
julia -t 2 --project=. scripts/gravity_lap.jl end1 64 > logs/end1.log 2>&1 &
P1=$!
julia -t 2 --project=. scripts/gravity_lap.jl end2 64 > logs/end2.log 2>&1 &
P2=$!
wait $P1; R1=$?
wait $P2; R2=$?
echo "end1 exit=$R1  end2 exit=$R2  at $(date +%T)"
[ $R1 -eq 0 ] && [ $R2 -eq 0 ] || { echo "a branch failed"; tail -5 logs/end1.log logs/end2.log; exit 1; }

julia --project=examples examples/plot_gravity_lap.jl 2>&1 | tail -8
