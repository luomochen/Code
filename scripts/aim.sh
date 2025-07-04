#!/bin/bash
#----------------------------------------
# Automiatically process AIM analysis.
#----------------------------------------
unzip CHGCAR.zip
chgsum.pl AECCAR0 AECCAR2
bader CHGCAR -ref CHGCAR_sum
#rm CHGCAR
echo "critic2 starts working."
cat > cps.cri <<!
crystal AECCAR0
load AECCAR0
load AECCAR2
load as "\$1+\$2" id rhoae
reference rhoae
MESHTYPE FRANCHINI VERYGOOD
AUTO SEED MESH 
cpreport cps.POSCAR
cpreport bond.cif graph
cpreport cps.json
!
critic2 cps.cri #> critic2.log
mv cps.POSCAR cps.vasp
if [ -f "ELFCAR" ]; then
    echo "File \"ELFCAR\" exists"
    mv ELFCAR ELFCAR.vasp
fi
#SEED SPHERE X0 0.99538 0.12738 0.08971 RADIUS 5 NPHI 10 NTHETA 10 NR 5 CLIP SPHERE 0.99538 0.12738 0.08971 5