#!/bin/bash
#-----------------------------------------------
# Use vaspkit to generate processed data of dos.
#-----------------------------------------------
ELEMENTS='Zr O'
ORBITALS=('s' 'py pz px' 'dxy dyz dz2 dxz dx2-y2')
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
echo 111 | vaspkit > /dev/null
echo TDOS.dat\&ITDOS.dat is generated.