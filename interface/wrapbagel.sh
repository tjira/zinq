#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -a <active> Number of active orbitals. (default: ${NACTIVE})
  -b <basis>  Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <closed> Number of closed orbitals. (default: ${NCLOSED})
  -g <charge> Charge of the system. (default: ${CHARGE})
  -k <nstate> Number of states to compute. (default: ${NSTATE})
  -m <method> Method to perform. (default: ${METHOD})
  -p <nspin>  The number associated with the spin states: 0 for singlet, 1 for doublet, 2 for triplet, etc. (default: ${NSPIN})
  -s <system> System file. (default: ${METHOD})
  -h          Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="svp"; CHARGE=0; NSPIN=0; METHOD="casscf"; NSTATE=2; NACTIVE=2; NCLOSED=4

while getopts "a:b:c:g:k:m:p:rs:h" OPT; do case "$OPT" in
  a ) NACTIVE="$OPTARG" ;;
  b ) BASIS="$OPTARG" ;;
  c ) NCLOSED="$OPTARG" ;;
  g ) CHARGE="$OPTARG" ;;
  k ) NSTATE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) NSPIN="$OPTARG" ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

cd $(dirname $SYSTEM) && SYSTEM=$(basename $SYSTEM)

METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM)

cat << EOT > bagel.json
{ "bagel" : [

{
  "title" : "molecule",
  "basis" : "$BASIS",
  "df_basis" : "$BASIS-jkfit",
  "angstrom" : true,
  "geometry" : [
EOT

awk 'NR == 1 {A = $1} NR > 2 {printf("{\"atom\" : \"%s\", \"xyz\" : [%20.14f, %20.14f, %20.14f]}%s\n", $1, $2, $3, $4, NR == A + 2 ? "" : ",")}' $SYSTEM >> bagel.json

cat << EOT >> bagel.json
  ]
},
EOT

[[ -f orbitals.archive ]] && cat << EOT >> bagel.json
{
  "title" : "load_ref",
  "file" : "orbitals",
  "continue_geom" : false
},
EOT

cat << EOT >> bagel.json
{
  "title" : "forces",
  "grads" : [
EOT

for I in $(seq 0 $((NSTATE - 1))); do
    awk -v I=$I -v NSTATE=$NSTATE 'BEGIN {printf("{\"title\" : \"force\", \"target\" : %d}%s", I, NSTATE > 1 ? "," : "")}' >> bagel.json
done

for I in $(seq 0 $((NSTATE - 1))); do
    for J in $(seq $((I + 1)) $((NSTATE - 1))); do
        awk -v I=$I -v J=$J -v NSTATE=$NSTATE 'BEGIN {printf("{\"title\" : \"nacme\", \"target\" : %d, \"target2\" : %d}%s", I, J, I == NSTATE - 2 && J == NSTATE - 1 ? "" : ",")}' >> bagel.json
    done
done

cat << EOT >> bagel.json
  ],
  "export" : true,
  "method" : [{
    "title" : "casscf",
    "charge" : $CHARGE,
    "nspin" : $NSPIN,
    "nact" : $NACTIVE,
    "nclosed" : $NCLOSED,
    "nstate" : $NSTATE,
    "maxiter" : 2000
  }]
},
{
  "title" : "save_ref",
  "file" : "orbitals"
}
]}
EOT

jq . bagel.json > bagel.json.tmp && mv bagel.json.tmp bagel.json && rm -f ENERGY.mat GRADIENT.mat NACV.mat

BAGEL bagel.json 2>&1 | tee bagel.out

if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print NSTATE, 1}'                                     >   ENERGY.mat
awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print 3 * NATOM * NSTATE, 1}'                         > GRADIENT.mat
awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print 3 * NATOM * (NSTATE * NSTATE - NSTATE) / 2, 1}' >     NACV.mat

awk '{printf("%20.14f\n", $1)}' ENERGY.out >> ENERGY.mat

for I in $(seq 0 $((NSTATE - 1))); do
    awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' FORCE_$I.out >> GRADIENT.mat
done

for I in $(seq 0 $((NSTATE - 1))); do
    for J in $(seq $((I + 1)) $((NSTATE - 1))); do
        awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' NACME_${I}_$J.out >> NACV.mat
    done
done

cat $SYSTEM >> $SYSTEM.all && cat bagel.out >> bagel.out.all
