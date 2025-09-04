#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "USAGE: $0 BASIS" >&2 &&  exit 1
fi

bse get-basis "$1" json --unc-spdf > "$1.json"
