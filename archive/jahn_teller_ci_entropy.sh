#!/bin/bash

D=10; MAX_CS=8; MAX_DE=0.2; DE_STEP=0.01; LEGEND=""

echo "$MAX_CS $((MAX_CID + 2))" > JAHN_TELLER_CI_ENTROPY.mat

for DE in $(seq 0 $DE_STEP $MAX_DE); do

    ROW=$(printf "%.2f" "$DE")

    for CS in $(seq 1 $MAX_CS); do

        CID=$(($D-$CS-1)); ITER=$(($CID*100000))

        cp jahn_teller_ci_entropy.json input.json

        jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.momentum_mean = [range($d) | 0   ]' input.json > tmp.json && mv tmp.json input.json
        jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.position_mean = [range($d) | 1e-5]' input.json > tmp.json && mv tmp.json input.json
        jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.gamma_mean    = [range($d) | 0   ]' input.json > tmp.json && mv tmp.json input.json
        jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.mass          = [range($d) | 1   ]' input.json > tmp.json && mv tmp.json input.json
        jq --argjson d "$D" '.zinq[0].options.potential.jahn_teller.d                = $d                ' input.json > tmp.json && mv tmp.json input.json
        jq --argjson d "$D" '.zinq[0].options.potential.jahn_teller.abs_coupling     = true              ' input.json > tmp.json && mv tmp.json input.json

        jq --argjson de "$DE" '.zinq[0].options.bias.variable.potential_energy_difference.value = $de' input.json > tmp.json && mv tmp.json input.json
        jq --argjson cs "$CS" '.zinq[0].options.potential.jahn_teller.coupling_sum              = $cs' input.json > tmp.json && mv tmp.json input.json

        jq --argjson cid "$CID" --argjson de "$DE" '.zinq[0].options.write.position_mean               = "POSITION_\($cid)D_DE=\($de).mat"' input.json > tmp.json && mv tmp.json input.json
        jq --argjson cid "$CID" --argjson de "$DE" '.zinq[0].options.write.state_potential_energy_mean = "PES_\($cid)D_DE=\($de).mat"     ' input.json > tmp.json && mv tmp.json input.json

        jq --argjson iter "$ITER" '.zinq[0].options.iterations = $iter' input.json > tmp.json && mv tmp.json input.json

        ENTROPY=$(zinq | grep ENTROPY | awk '{print $6}')

        ROW=$(printf "%s %12.8f" "$ROW" "$ENTROPY")

        printf "%03d %.2f %12.8f\n" "$CID" "$DE" "$ENTROPY"
    done

    echo "$ROW" >> JAHN_TELLER_CI_ENTROPY.mat
done

rm input.json

plot JAHN_TELLER_CI_ENTROPY.mat --legends every "8D Seam" "7D Seam" "6D Seam" "5D Seam" "4D Seam" "3D Seam" "2D Seam" "1D Seam" --title "Seam Entropy as a Function of $\Delta E$ on a ${D}D Potential" --xlabel "$\Delta E$ (Hartree)" --ylabel "Entropy of CI (J/K/Mol)" --output JAHN_TELLER_CI.png
