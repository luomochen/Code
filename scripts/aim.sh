#!/bin/bash
#----------------------------------------
# Automiatically process AIM analysis.
#----------------------------------------
unzip CHGCAR.zip
chgsum.pl AECCAR0 AECCAR2
bader CHGCAR -ref CHGCAR_sum
rm CHGCAR
cat > cps.cri <<!
crystal AECCAR0
load AECCAR0
load AECCAR2
load as "\$1+\$2" id rhoae
reference rhoae
auto seed pair seed ws seed triplet
cpreport cps.POSCAR
cpreport bond.cif graph
cpreport cps.json
!
critic2 cps.cri > critic2.log
mv cps.POSCAR cps.vasp
if [ -f "ELFCAR" ]; then
    echo "File \"ELFCAR\" exists"
    mv ELFCAR ELFCAR.vasp
fi
