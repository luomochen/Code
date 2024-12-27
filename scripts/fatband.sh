#!/bin/bash
#---------------------------------------------------
# Use vaspkit to generate processed data of fatband.
#---------------------------------------------------
declare -A ORBITALS=(
    [1]="s"
    [2]="py"
    [3]="pz"
    [4]="px"
    [5]="dxy"
    [6]="dyz"
    [7]="dz2"
    [8]="dxz"
    [9]="dx2-y2"
)

ELEMENTS="Zr O"
for E in $ELEMENTS
do
    for i in {1..4}; do
        echo -e "216\n1\n$E\n$i\n" | vaspkit > /dev/null
        if [ -f "PBAND_SUM.dat" ]; then
            mv PBAND_SUM.dat "$E"_"${ORBITALS[$i]}".dat
            echo "$E"_"${ORBITALS[$i]}".dat is generated.
        fi 
    done 
done