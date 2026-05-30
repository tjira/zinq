#!/bin/bash

perf record -g -- "$@"

perf script | gprof2dot -e 1 -f perf -n 5 -z "main.main" | dot -T svg -o profile.svg

rm -f perf.data
