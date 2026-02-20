#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -b <basis>        Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <charge>       Charge of the system. (default: ${CHARGE})
  -m <method>       Method to perform. (default: ${METHOD})
  -p <multiplicity> Spin multiplicity of the system. (default: ${MULTIPLICITY})
  -s <system>       System file. (default: ${SYSTEM})
  -h                Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="sto-3g"; CHARGE=0; MULTIPLICITY=1; METHOD="hf";

while getopts "b:c:k:m:p:rs:h" OPT; do case "$OPT" in
  b ) BASIS="$OPTARG" ;;
  c ) CHARGE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) MULTIPLICITY="$OPTARG" ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

cd $(dirname $SYSTEM) && SYSTEM=$(basename $SYSTEM)

METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM); GUESS=$([ -f orca.gbw ] && echo "MOREAD" || echo "HCORE")

if [[ $METHOD == "FCI" ]]; then
    METHOD=""; OPTIONS="%CASSCF DOFCI TRUE; END"
fi

if [[ -f orca.gbw ]]; then
    mv orca.gbw guess.gbw
fi

cat << EOT > orca.inp
! $METHOD ${BASIS^^} $GUESS KDIIS NOFROZENCORE TIGHTSCF ENGRAD LARGEPRINT

%moinp "guess.gbw"

*xyzfile $CHARGE $MULTIPLICITY $SYSTEM
EOT

sed -i 's/   /  /g ; s/  / /g' orca.inp && rm -f ENERGY.mat GRADIENT.mat NACV.mat

orca orca.inp 2>&1 | tee orca.out

if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

awk -v NATOM=$NATOM 'BEGIN {print 1, 1} NR == 8 {printf "%20.14f\n", $1}' orca.engrad > ENERGY.mat

awk -v NATOM=$NATOM 'BEGIN {print 3 * NATOM, 1} NR > 11 && NR <= 3 * NATOM + 11 {printf "%20.14f\n", $1}' orca.engrad > GRADIENT.mat
