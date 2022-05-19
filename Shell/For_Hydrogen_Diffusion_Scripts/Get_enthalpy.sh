#!/bin/bash

# 该脚本用于获取焓值。
# 将文件夹名称存储在一个变量中。
list="P21c Pbca Pnma P62m P63mmc"
# 分别进入不同的文件夹中。
for spacegroup in $list
do 
    cd /d/"OneDrive - lzu.edu.cn"/Hydrogen_Diffusion_Oxides/relax_test3||exit
    cd "$spacegroup"||exit
    # 由于0 GPa下PV = 0，所以内能等于焓值。
    # 利用管道过滤出最后一行值后重定向至输出文件。
    grep without "$spacegroup"_0GPa/OUTCAR | tail -n 1 | gawk '{print $4}' >> "$spacegroup"_enthalpy
    for ((i = 1; i < 29; i++))
    do
        if [ "$i" -lt 11 ]
        then
            a=$[$i * 5]
            grep enthalpy "$spacegroup"_"$a"GPa/OUTCAR | tail -n 1 | gawk '{print $5}' >> "$spacegroup"_enthalpy 
        # 加压间隔为25GPa，逐步加压至500GPa。
        else
            if [ "$i" = 11 ]
            then 
                grep enthalpy "$spacegroup"_75GPa/OUTCAR | tail -n 1 | gawk '{print $5}' >> "$spacegroup"_enthalpy 
            else
                a=$[($i - 10) * 25 + 50]
                grep enthalpy "$spacegroup"_"$a"GPa/OUTCAR | tail -n 1 | gawk '{print $5}' >> "$spacegroup"_enthalpy
            fi
        fi
    done
done 
