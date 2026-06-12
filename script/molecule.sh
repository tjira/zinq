#!/bin/bash

rm -rf temp && mkdir -p example/molecule && BASIS="STO-3G"; METHOD="RHF"

mkdir -p temp && cat > temp/CO2.xyz << EOM
3
CO2
C  0.0000 0.0000 0.0000
O -1.0000 0.0000 0.0000
O  1.0000 0.0000 0.0000
EOM

mkdir -p temp && cat > temp/H2.xyz << EOM
2
H2
H 0.0000 0.0000 0.0000
H 1.0000 0.0000 0.0000
EOM

mkdir -p temp && cat > temp/HCl.xyz << EOM
2
HCl
H  0.0000 0.0000 0.0000
Cl 1.0000 0.0000 0.0000
EOM

mkdir -p temp && cat > temp/HF.xyz << EOM
2
HF
H 0.0000 0.0000 0.0000
F 1.0000 0.0000 0.0000
EOM

mkdir -p temp && cat > temp/Ne2.xyz << EOM
2
Ne2
Ne 0.0000 0.0000 0.0000
Ne 1.0000 0.0000 0.0000
EOM

mkdir -p temp && cat > temp/acetone.xyz << EOM
10
acetone
O -0.4797 -0.0446  1.2294
C  1.4376  0.0818 -0.1875
C -0.0377 -0.0035  0.0968
C -0.9294 -0.0345 -1.1150
H  2.0471  0.1029  0.7251
H  1.6643  0.9894 -0.7620
H  1.7616 -0.7796 -0.7860
H -0.7844  0.8703 -1.7195
H -0.6856 -0.8986 -1.7467
H -1.9948 -0.0979 -0.8587
EOM

mkdir -p temp && cat > temp/water.xyz << EOM
3
water
O -0.0669 -0.0232 -0.0267
H  0.6848 -0.6120  0.5191
H  0.3783  0.9803 -0.0938
EOM

# PROCESS EACH XYZ FILE
cd temp && for MOL in *.xyz; do

    # CREATE ORCA INPUT FILE
    echo -e "! $METHOD $BASIS OPT\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"

    # EXTRACT ENERGY
    ENERGY=$(grep "FINAL SINGLE POINT ENERGY" "${MOL%.*}.out" | tail -n 1 | awk '{print $5}')

    # UPDATE XYZ FILE
    sed -i "2s/.*/${MOL%.*}\/$METHOD\/$BASIS\/$ENERGY/" "$MOL"

    # LOOP THROUGH XYZ FILE
    while read -r LINE; do

        # FORMAT THE FIRST TWO LINES
        echo "$LINE" | awk 'NF==1 {printf("%s\n", $1)}' >> "$MOL.tmp"

        # FORMAT LINES WITH COORDINATES
        echo "$LINE" | awk 'NF==4 {printf("%2s % 3.8f % 3.8f % 3.8f\n", $1, $2, $3, $4)}' >> "$MOL.tmp"

    # PIPE THE XYZ FILE TO THE WHILE LOOP
    done < "$MOL"

    # MOVE THE TEMP FILE TO THE EXAMPLE DIRECTORY
    mv "$MOL.tmp" "$MOL" && mv "$MOL" "../example/molecule/$MOL"

# FINISH THE FOR LOOP OVER ALL XYZ FILES
done && cd .. && rm -rf temp
