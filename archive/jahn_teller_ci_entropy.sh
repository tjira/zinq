#!/bin/bash

MAX_D=100

echo "$MAX_D 2" > JAHN_TELLER_CI_ENTROPY.mat

for D in $(seq 2 $MAX_D); do

    cp jahn_teller_ci_entropy.json input.json

    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.momentum_mean = [range($d) | 0   ]   ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.position_mean = [range($d) | 1e-5]   ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.gamma_mean    = [range($d) | 0   ]   ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.initial_conditions.model.mass          = [range($d) | 1   ]   ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.potential.jahn_teller.d                = $d                   ' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.write.position_mean                    = "POSITION_\($d)D.mat"' input.json > tmp.json && mv tmp.json input.json
    jq --argjson d "$D" '.zinq[0].options.write.state_potential_energy_mean      = "PES_\($d)D.mat"     ' input.json > tmp.json && mv tmp.json input.json

    ENTROPY=$(zinq | grep SCHLITTER | awk '{print $6}')

    printf "%03d %12.8f\n" "$D" "$ENTROPY" >> JAHN_TELLER_CI_ENTROPY.mat

    tail -n 1 JAHN_TELLER_CI_ENTROPY.mat
done

rm input.json

plot JAHN_TELLER_CI_ENTROPY.mat POSITION_4D.mat:2,3 --figsize 6 16 --subplots 121 122 122 122 122 --legends 1,2 "First CI Coordinate" "Second CI Coordinate" --title "Schlitter Entropy of CI as a Function of Dimension" "CI Coordinates of 4D Model as a function of Time" --xlabel "Dimension" "Time (a.u.)" --ylabel "Schlitter Entropy of CI (J/K/Mol)" "Position (a.u.)" --output JAHN_TELLER_CI.png
