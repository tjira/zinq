#!/bin/bash

perf record -g -- "$@"

perf script | gprof2dot -f perf | dot -T pdf -o profile.pdf

rm -f perf.data
