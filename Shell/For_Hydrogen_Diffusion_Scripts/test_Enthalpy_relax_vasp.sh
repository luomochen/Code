#!/bin/bash

# 该脚本用于对Enthlapy_relax_vasp.pbs进行测试。
# 将文件夹名称存储在一个变量中。
list="P21c Pbca Pnma P62m P63mmc"
# 分别进入不同的文件夹中。
for spacegroup in $list
do  
    # 优化未加压的结构，获得CONTCAR后再进入循环。
    # 进入存放输入文件的relax文件夹。
    cd /path||exit
    cd "$spacegroup"||exit
    cd "$spacegroup"_0GPa||exit
    mpirun -np 4 vasp
    cd ..
    for ((i = 1; i < 24; i++))
    do
        # 加压间隔为10GPa，逐步加至50GPa。
        if [ "$i" -lt 6 ]
        then
            a=$[$i * 1000000]
            b=$[($i - 1) * 10]
            c=$[$i * 10]

            mkdir "$spacegroup"_"$c"GPa
            cp "$spacegroup"_"$b"GPa/INCAR "$spacegroup"_"$b"GPa/POTCAR "$spacegroup"_"$b"GPa/KPOINTS "$spacegroup"_"$c"GPa
            cp "$spacegroup"_"$b"GPa/CONTCAR "$spacegroup"_"$c"GPa/POSCAR
            sed -i "20c PSTRESS = $a" "$spacegroup"_"$c"GPa/INCAR
            cd "$spacegroup"_"$c"GPa||exit
            mpirun -np 4 vasp
            cd ..
        # 加压间隔为25GPa，逐步加压至500GPa。
        else
            if [ "$i" = 6 ]
            then 
                mkdir "$spacegroup"_75GPa
                cp "$spacegroup"_50GPa/INCAR "$spacegroup"_50GPa/POTCAR "$spacegroup"_50GPa/KPOINTS "$spacegroup"_75GPa
                cp "$spacegroup"_50GPa/CONTCAR "$spacegroup"_75GPa/POSCAR
                sed -i "20c PSTRESS = 7500000" "$spacegroup"_75GPa/INCAR
                cd "$spacegroup"_75GPa||exit
                mpirun -np 4 vasp
                cd ..
            else
                a=$[($i - 5) * 2500000 + 5000000]
                b=$[($i - 6) * 25 + 50]
                c=$[($i - 5) * 25 + 50]

                mkdir "$spacegroup"_"$c"GPa
                cp "$spacegroup"_"$b"GPa/INCAR "$spacegroup"_"$b"GPa/POTCAR "$spacegroup"_"$b"GPa/KPOINTS "$spacegroup"_"$c"GPa
                cp "$spacegroup"_"$b"GPa/CONTCAR "$spacegroup"_"$c"GPa/POSCAR
                sed -i "20c PSTRESS = $a" "$spacegroup"_"$c"GPa/INCAR
                cd "$spacegroup"_"$c"GPa||exit
                mpirun -np 4 vasp
                cd ..
            fi
        fi
    done
done