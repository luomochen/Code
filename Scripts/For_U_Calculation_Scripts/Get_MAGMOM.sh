#!/bin/bash

#touch MAGMOM
for((i=0;i<32;i++)); do
    k=$(echo $i+9 | bc)
    a=$(sed -n "$k"p POSCAR.vasp | gawk '{print $3}' | tr -d $'\r')
    D1=$(echo "$a - 0.750000000" | bc)
    D1=${D1#-}
    D1=$(echo "$D1 < 0.1" | bc)
    D2=$(echo "$a - 0.250000000" | bc)
    D2=${D2#-}
    D2=$(echo "$D2 < 0.1" | bc)
    if [ "$D1" -eq 1 ] || [ "$D2" -eq 1 ]
    then
        echo -n "2.6 " >> MAGMOM
    else
        echo -n "-2.6 " >> MAGMOM
    fi
done