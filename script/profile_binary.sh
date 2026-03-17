#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: $0 <COMMAND>" && exit 1
fi

valgrind --callgrind-out-file=callgrind.out --tool=callgrind $1 && gprof2dot -e 1 -f callgrind -n 5 -o profile.dot callgrind.out && dot -T pdf profile.dot -o profile.pdf
