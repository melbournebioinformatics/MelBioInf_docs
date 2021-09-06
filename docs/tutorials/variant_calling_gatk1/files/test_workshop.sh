#!/bin/bash

N=5
seq -w 1 "$N" | parallel --joblog test_workshop.log -j 0 --halt 2 'RUNNUMBER={} ./run_complete_workshop.sh'
