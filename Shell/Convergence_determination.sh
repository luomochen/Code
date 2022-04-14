#!/bin/bash

# 该脚本用于确定批量优化的结构是否收敛。
# 将文件夹名称存储在一个变量中。
list="P21c Pbca Pnma P62m P63mmc"
# 分别进入不同的文件夹中。
for spacegroup in $list
do  
    cd /d/"OneDrive - lzu.edu.cn"/Hydro_Diff_Oxi/Experiment/relax_test10||exit
    cd "$spacegroup"||exit
    touch "$spacegroup"_convergence_info
    for ((i = 0; i < 24; i++))
    do
        if [ "$i" -lt 6 ]
        then
            c=$[$i * 10]

            cd "$spacegroup"_"$c"GPa||exit
            if grep "required accuracy" OUTCAR
            then
                echo "$c Convergence" >> ../"$spacegroup"_convergence_info
            else
                echo "$c Non-Convergence" >> ../"$spacegroup"_convergence_info
            fi
            cd ..
        else
            if [ "$i" = 6 ]
            then 
                cd "$spacegroup"_75GPa||exit
                if grep "required accuracy" OUTCAR
                then
                    echo "$c Convergence" >> ../"$spacegroup"_convergence_info
                else
                    echo "$c Non-Convergence" >> ../"$spacegroup"_convergence_info
                fi
                cd ..
            else
                c=$[($i - 5) * 25 + 50]

                cd "$spacegroup"_"$c"GPa||exit
                if grep "required accuracy" OUTCAR
                then
                    echo "$c Convergence" >> ../"$spacegroup"_convergence_info
                else
                    echo "$c Non-Convergence" >> ../"$spacegroup"_convergence_info
                fi
                cd ..
            fi
        fi
    done
done
