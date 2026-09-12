#!/bin/zsh
# Sequential long runs (each ~15–60 min). Logs in logs/, CSVs in out/.
cd ~/Developer/FallingPlate2D.jl
run() { name=$1; shift; echo "=== $(date '+%H:%M') start $name" >> logs/chain.log
        julia -t 4 --project=. examples/free_fall.jl "$@" $name > logs/run_$name.log 2>&1
        echo "=== $(date '+%H:%M') end $name exit $?" >> logs/chain.log; }
run andersen_dns_b8_L96   8 2.52 2090 96  60
run andersen_Da1e-4_L64   8 2.52 2090 64  60 1e-4
run andersen_Da1e-3_L64   8 2.52 2090 64  60 1e-3
run andersen_Da1e-2_L64   8 2.52 2090 64  60 1e-2
run andersen_dns_b8_L128  8 2.52 2090 128 60
echo "=== $(date '+%H:%M') CHAIN DONE" >> logs/chain.log
