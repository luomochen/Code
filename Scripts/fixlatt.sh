#!/bin/bash
#---------------------------------------------------------------
# This script is used to ensure that the lattice constants of 
# the structure file do not change.
#---------------------------------------------------------------
if [ -f "../pft.vasp" ]; then
    a=$(sed -n "3p" ../pft.vasp)
    b=$(sed -n "4p" ../pft.vasp)
    c=$(sed -n "5p" ../pft.vasp)
    list=$( ls -d */ )
    for dir in $list; do
        if [[ $dir == 0* ]] 
        then
            if [ -f "$dir/POSCAR" ]; then
                sed -i "3c \ $a" "$dir/POSCAR"
                sed -i "4c \ $b" "$dir/POSCAR"
                sed -i "5c \ $c" "$dir/POSCAR"
            fi
            if [ -f "$dir/CONTCAR" ]; then
                sed -i "3c \ $a" "$dir/CONTCAR"
                sed -i "4c \ $b" "$dir/CONTCAR"
                sed -i "5c \ $c" "$dir/CONTCAR"
            fi          
        fi
    done
    read -rp "Do you want to change the lattice constants of the states?[y/n]: " A
    if [ "$A" == y ]; then
        list="ini fin"
        for dir in $list; do
            sed -i "3c \ $a" "$dir/POSCAR"
            sed -i "4c \ $b" "$dir/POSCAR"
            sed -i "5c \ $c" "$dir/POSCAR"            
            sed -i "3c \ $a" "$dir/CONTCAR"
            sed -i "4c \ $b" "$dir/CONTCAR"
            sed -i "5c \ $c" "$dir/CONTCAR"
        done
    fi
    echo "Lattice constants changed!"    
else
   echo "Perfect crystal structure file not found!" 
fi 