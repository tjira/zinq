#!/bin/bash

perf record -g -- "$@"

perf script | gprof2dot -e 1 -f perf -n 5 | dot -T pdf -o profile.pdf

rm -f perf.data
