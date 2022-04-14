#!/bin/bash

# 获取优化后的氢的点位和系统的总能量。
list="m o1 o2"
for phase in $list
do
    cd /d/"OneDrive - lzu.edu.cn"/Hydro_Diff_Oxi/Experiment/H_optim_sites/test2||exit
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
                cd "$phase"_"$i""$j""$k"||exit
                r1=$(sed -n '7p' CONTCAR | gawk '{print $1}')
                r2=$(sed -n '7p' CONTCAR | gawk '{print $2}')
                r3=$(sed -n '7p' CONTCAR | gawk '{print $3}')
                r=$["$r1" + "$r2" + "$r3" + 9]
                coordinate=$(sed -n "$r"p CONTCAR | gawk '{print $1" "$2" "$3}')
                energy=$(grep without OUTCAR | tail -n 1 | gawk '{print $4}')
                echo "$a $b $c $coordinate $energy" >> "$phase"_data
                cd ..
            done
        done
    done
done
