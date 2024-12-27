#!/bin/bash
#-----------------------------------------------
# Use vaspkit to generate processed data of dos.
#-----------------------------------------------
ELEMENTS='Zr O'
ORBITALS=('s' 'py pz px')
for E in $ELEMENTS
do
    for O in "${ORBITALS[@]}"; do
        echo -e "115\n$E\n$O\n" | vaspkit > /dev/null
        if [ -f "PDOS_USER.dat" ]; then
            O=${O// /_}
            mv PDOS_USER.dat "$E"_"$O".dat
            echo "$E"_"$O".dat is generated.
        fi
    done 
done