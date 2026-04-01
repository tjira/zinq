#!/bin/bash

D=20; MAX_CID=19

echo "$MAX_CID 2" > JAHN_TELLER_CI_ENTROPY.mat

for CS in $(seq 0 $MAX_CID); do

    CID=$(($D-$CS-1)); cp jahn_teller_ci_entropy.json input.json

    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.momentum_mean = [range($d) | 0   ]' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.position_mean = [range($d) | 1e-5]' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.gamma_mean    = [range($d) | 0   ]' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.mass          = [range($d) | 1   ]' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.potential.jahn_teller.d                = $d                ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.potential.jahn_teller.abs_coupling     = true              ' input.json > tmp.json && mv tmp.json input.json

    jq --argjson cs "$CS" '.zinq[0].options.potential.jahn_teller.coupling_sum = $cs' input.json > tmp.json && mv tmp.json input.json

    jq --argjson cid "$CID" '.zinq[0].options.write.position_mean                = "POSITION_\($cid)D.mat"' input.json > tmp.json && mv tmp.json input.json
    jq --argjson cid "$CID" '.zinq[0].options.write.state_potential_energy_mean  = "PES_\($cid)D.mat"     ' input.json > tmp.json && mv tmp.json input.json

    ENTROPY=$(zinq | grep SCHLITTER | awk '{print $6}')

    printf "%03d %12.8f\n" "$CID" "$ENTROPY" >> JAHN_TELLER_CI_ENTROPY.mat

    tail -n 1 JAHN_TELLER_CI_ENTROPY.mat
done

rm input.json

plot JAHN_TELLER_CI_ENTROPY.mat --title "Seam Entropy as a Function of Dimension on a ${D}D Potential" --xlabel "Seam Dimension" --ylabel "Schlitter Entropy of CI (J/K/Mol)" --output JAHN_TELLER_CI.png
