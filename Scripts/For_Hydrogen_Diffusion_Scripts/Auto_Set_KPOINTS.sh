#!/bin/bash
#---------------------------------------------------------------
# 该脚本用于随晶格常数变化自动选取合适的N_1,N_2和N_3并生成KPOINTS。
#---------------------------------------------------------------

# 生成KPOINTS文件头部。
cat > KPOINTS <<!
KPOINTS
 0
Gamma
!

# 根据经验Ka=Kb=Kc=30~40。
a1=$(sed -n '3p' POSCAR | gawk '{print $1}')
a2=$(sed -n '3p' POSCAR | gawk '{print $2}')
a3=$(sed -n '3p' POSCAR | gawk '{print $3}')
a=$(echo "sqrt($a1^2 + $a2^2 + $a3^2)" | bc)
N_1=$(echo "scale = 0; 40 / $a" | bc)

b1=$(sed -n '4p' POSCAR | gawk '{print $1}')
b2=$(sed -n '4p' POSCAR | gawk '{print $2}')
b3=$(sed -n '4p' POSCAR | gawk '{print $3}')
b=$(echo "sqrt($b1^2 + $b2^2 + $b3^2)" | bc)
N_2=$(echo "scale = 0; 40 / $b" | bc)

c1=$(sed -n '5p' POSCAR | gawk '{print $1}')
c2=$(sed -n '5p' POSCAR | gawk '{print $2}')
c3=$(sed -n '5p' POSCAR | gawk '{print $3}')
c=$(echo "sqrt($c1^2 + $c2^2 + $c3^2)" | bc)
N_3=$(echo "scale = 0; 40 / $c" | bc)

cat >> KPOINTS <<!
$N_1 $N_2 $N_3
0 0 0
!