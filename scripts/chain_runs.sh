#!/bin/zsh
# Sequential long runs (each ~15–60 min). Logs in logs/, CSVs in out/.
# usage: run NAME χ ρ Re L tend [Da]
cd ~/Developer/FallingPlate2D.jl
run() { name=$1; shift; echo "=== $(date '+%H:%M') start $name" >> logs/chain.log
        julia -t 4 --project=. examples/free_fall.jl "$1" "$2" "$3" "$4" "$5" "$name" "${6:-0}" > logs/run_$name.log 2>&1
        echo "=== $(date '+%H:%M') end $name exit $?" >> logs/chain.log; }
if [[ "$1" == "after" ]]; then   # wait for a previous chain to finish first
    until grep -q "CHAIN DONE" logs/chain.log; do sleep 30; done
fi
run andersen_Da1e-4_L64   8 2.52 2090 64 60 1e-4
run andersen_Da1e-3_L64   8 2.52 2090 64 60 1e-3
run andersen_Da1e-2_L64   8 2.52 2090 64 60 1e-2
run andersen_Re1500_L64   8 2.52 1500 64 60
echo "=== $(date '+%H:%M') CHAIN2 DONE" >> logs/chain.log
