#!/bin/zsh
# Do NOT edit while running (zsh reads scripts incrementally).
cd ~/Developer/FallingPlate2D.jl
echo "=== $(date '+%H:%M') start L=64 set" >> logs/chain_gust.log
julia -t 4 --project=. scripts/gust_runs.jl 64 > logs/gust_L64.log 2>&1
echo "=== $(date '+%H:%M') end L=64 exit $?" >> logs/chain_gust.log
echo "=== $(date '+%H:%M') start L=96 set" >> logs/chain_gust.log
julia -t 4 --project=. scripts/gust_runs.jl 96 > logs/gust_L96.log 2>&1
echo "=== $(date '+%H:%M') end L=96 exit $?" >> logs/chain_gust.log
echo "=== $(date '+%H:%M') GUST CHAIN DONE" >> logs/chain_gust.log
