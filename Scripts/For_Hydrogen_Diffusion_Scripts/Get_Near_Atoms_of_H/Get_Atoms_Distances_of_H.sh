#!/bin/bash

# 该脚本用于确定过渡态计算中距离扩散原子距离较近的原子，确定结构的合理程度。

mkdir Atom_distances
for((i=1;i<6;i++));
do
    echo "Atom Number Distance" >> Atom_distances/Atom_distances_0$i.txt
    # 读取晶胞基矢，用于之后的分数坐标向笛卡尔坐标转换。
    a1=$(sed -n '3p' 0$i/CONTCAR | gawk '{print $1}')
    a2=$(sed -n '3p' 0$i/CONTCAR | gawk '{print $2}')
    a3=$(sed -n '3p' 0$i/CONTCAR | gawk '{print $3}')

    b1=$(sed -n '4p' 0$i/CONTCAR | gawk '{print $1}')
    b2=$(sed -n '4p' 0$i/CONTCAR | gawk '{print $2}')
    b3=$(sed -n '4p' 0$i/CONTCAR | gawk '{print $3}')

    c1=$(sed -n '5p' 0$i/CONTCAR | gawk '{print $1}')
    c2=$(sed -n '5p' 0$i/CONTCAR | gawk '{print $2}')
    c3=$(sed -n '5p' 0$i/CONTCAR | gawk '{print $3}')

    # 读取氢原子坐标。
    H_x_coordinate=$(sed -n '105p' 0$i/CONTCAR | gawk '{print $1}')
    H_y_coordinate=$(sed -n '105p' 0$i/CONTCAR | gawk '{print $2}')
    H_z_coordinate=$(sed -n '105p' 0$i/CONTCAR | gawk '{print $3}')

    # 确定原子种类和原子标号。
    Zr_number=$(sed -n '7p' 0$i/CONTCAR | gawk '{print $1}')
    O_number=$(sed -n '7p' 0$i/CONTCAR | gawk '{print $2}')
    Cell_atom_number=$Zr_number+$O_number+1

    # 开始读取每个CONTCAR中的原子坐标。
    for((j=0;j<Cell_atom_number;j++));
    do
        line_number=$(echo "$j+9" | bc)
        if [ "$(echo "$j+1 <= $Zr_number"|bc)" -eq 1 ] 
        then
            Element="Zr"
            serial_number=$(echo "$j+1" | bc)
        elif [ "$(echo "$j < $Cell_atom_number-1" |bc)" -eq 1 ] 
        then
            Element="O"
            serial_number=$(echo "$j+1-$Zr_number" | bc)
        else
            Element="H"
            serial_number=1
        fi
        x_coordinate=$(sed -n "$line_number"p 0$i/CONTCAR | gawk '{print $1}')
        y_coordinate=$(sed -n "$line_number"p 0$i/CONTCAR | gawk '{print $2}')
        z_coordinate=$(sed -n "$line_number"p 0$i/CONTCAR | gawk '{print $3}')
        # 开始计算该过渡态下每个Zr，O原子和氢原子之间的距离。
        x_direct_distance=$(echo "$x_coordinate - ${H_x_coordinate[$k]}" | bc)
        y_direct_distance=$(echo "$y_coordinate - ${H_y_coordinate[$k]}" | bc)
        z_direct_distance=$(echo "$z_coordinate - ${H_z_coordinate[$k]}" | bc)

        a_square=$(echo "($x_direct_distance*$a1 + $y_direct_distance*$b1 + $z_direct_distance*$c1)^2" | bc)
        b_square=$(echo "($x_direct_distance*$a2 + $y_direct_distance*$b2 + $z_direct_distance*$c2)^2" | bc)
        c_square=$(echo "($x_direct_distance*$a3 + $y_direct_distance*$b3 + $z_direct_distance*$c3)^2" | bc)

        distance=$(echo "sqrt($a_square + $b_square + $c_square)" | bc)
        echo "$Element $serial_number $distance" >> Atom_distances/Atom_distances_0$i.txt
    done
    # 向python脚本传递image编号参数，并调用python脚本进行排序。
    python Sort_Atoms_with_Distance.py $i
done