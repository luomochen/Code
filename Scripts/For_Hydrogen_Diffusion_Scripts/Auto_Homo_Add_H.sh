#!/bin/bash

# 该脚本用于对Auto_homo_add_Hydrogen.pbs进行测试。
# 该脚本用于在分数坐标下均匀地将单个氢原子均匀地分散在晶胞各处。
list="m o1 o2"
for phase in $list
do    
    cd /path||exit
    cd "$phase"||exit
    for ((i = 0; i < 5; i++))
    do 
        a=$(echo "scale = 17 ; $i * 0.2" | bc | awk '{printf "%.17f", $0}')
        for ((j = 0; j < 5; j++))
        do
            b=$(echo "scale = 17 ; $j * 0.2" | bc | awk '{printf "%.17f", $0}')
            for ((k = 0; k < 5; k++))
            do  
                c=$(echo "scale = 17 ; $k * 0.2" | bc | awk '{printf "%.17f", $0}')
                cp -r "$phase"_origin "$phase"_"$i""$j""$k"
                cd "$phase"_"$i""$j""$k"||exit
                echo " $a $b $c T T T" >> POSCAR
                mpirun -np 4 vasp                
                cd ..
            done
        done
    done
done
