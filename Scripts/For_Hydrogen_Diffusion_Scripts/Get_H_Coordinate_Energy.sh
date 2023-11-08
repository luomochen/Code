#!/bin/bash
#-------------------------------------------------------
# 获取优化后的氢的坐标和系统的总能量。
#-------------------------------------------------------

echo "x,y,z,energy" > Outcome.txt
echo "AA,AA,AA,eV" >> Outcome.txt
for ((i = 0; i < 5; i++))
do 
for ((j = 0; j < 5; j++))
do
for ((k = 0; k < 5; k++))
do  
cd "$i""$j""$k"||exit
r1=$(sed -n '7p' CONTCAR | gawk '{print $1}')
r2=$(sed -n '7p' CONTCAR | gawk '{print $2}')
r3=$(sed -n '7p' CONTCAR | gawk '{print $3}')
r=$(("$r1" + "$r2" + "$r3" + 9))
x=$(sed -n "$r"p CONTCAR | gawk '{print $1}')
y=$(sed -n "$r"p CONTCAR | gawk '{print $2}')
z=$(sed -n "$r"p CONTCAR | gawk '{print $3}')
energy=$(grep without OUTCAR | tail -n 1 | gawk '{print $4}')
echo "$x,$y,$z,$energy" >> Outcome.txt
cd ..
done
done
done