#!/bin/bash
#------------------------------------------------------------------
# 该脚本用于将远离扩散原子的原子固定，便于快速搜索过渡态。
#------------------------------------------------------------------

# 读取晶胞基矢，用于之后的分数坐标向笛卡尔坐标转换。
a1=$(sed -n '3p' 00/POSCAR | gawk '{print $1}')
a2=$(sed -n '3p' 00/POSCAR | gawk '{print $2}')
a3=$(sed -n '3p' 00/POSCAR | gawk '{print $3}')

b1=$(sed -n '4p' 00/POSCAR | gawk '{print $1}')
b2=$(sed -n '4p' 00/POSCAR | gawk '{print $2}')
b3=$(sed -n '4p' 00/POSCAR | gawk '{print $3}')

c1=$(sed -n '5p' 00/POSCAR | gawk '{print $1}')
c2=$(sed -n '5p' 00/POSCAR | gawk '{print $2}')
c3=$(sed -n '5p' 00/POSCAR | gawk '{print $3}')
# 设置空的数组用来存放每个过渡态中氢原子的坐标。
H_x_coordinate=()
H_y_coordinate=()
H_z_coordinate=()
for((i=1;i<6;i++));
do
    cp 0$i/POSCAR 0$i/POSCAR_unfix
    sed -i 's/ *$//g' 0$i/POSCAR
    sed -i "7a\Selective dynamics" 0$i/POSCAR
    H_x_coordinate[i]=$(sed -n '$p' 0$i/POSCAR | gawk '{print $1}')
    H_y_coordinate[i]=$(sed -n '$p' 0$i/POSCAR | gawk '{print $2}')
    H_z_coordinate[i]=$(sed -n '$p' 0$i/POSCAR | gawk '{print $3}')
done

Zr_number=$(sed -n '7p' 00/POSCAR | gawk '{print $1}')
O_number=$(sed -n '7p' 00/POSCAR | gawk '{print $2}')
Cell_atom_number=$Zr_number+$O_number+1
# 进入每个过渡态文件夹中。
for((i=1;i<6;i++));
do
    cd 0$i || exit
    # 开始读取每个POSCAR中的原子坐标。
    for((j=0;j<Cell_atom_number;j++));
    do
    line_number=$(echo "$j+10" | bc)
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
    x_coordinate=$(sed -n "$line_number"p POSCAR | gawk '{print $1}')
    y_coordinate=$(sed -n "$line_number"p POSCAR | gawk '{print $2}')
    z_coordinate=$(sed -n "$line_number"p POSCAR | gawk '{print $3}')
    # 开始计算该过渡态下每个Zr，O原子和所有过渡态下氢原子之间的距离，以确定扩散通道上所有可动原子。
    for((k=1;k<6;k++));
    do
        x_direct_distance=$(echo "$x_coordinate - ${H_x_coordinate[$k]}" | bc)
        y_direct_distance=$(echo "$y_coordinate - ${H_y_coordinate[$k]}" | bc)
        z_direct_distance=$(echo "$z_coordinate - ${H_z_coordinate[$k]}" | bc)

        a_square=$(echo "($x_direct_distance*$a1 + $y_direct_distance*$b1 + $z_direct_distance*$c1)^2" | bc)
        b_square=$(echo "($x_direct_distance*$a2 + $y_direct_distance*$b2 + $z_direct_distance*$c2)^2" | bc)
        c_square=$(echo "($x_direct_distance*$a3 + $y_direct_distance*$b3 + $z_direct_distance*$c3)^2" | bc)

        distance=$(echo "sqrt($a_square + $b_square + $c_square)" | bc)            
        # 只要H原子和周围原子的距离小于2埃，立刻固定并跳出循环，否则继续寻找。
        if [ "$(echo "$distance <= 2.0"|bc)" -eq 1 ] 
        then
            sed -i "$line_number"'s/$/& T T T/' POSCAR
            echo "$i $Element $serial_number $x_coordinate $y_coordinate $z_coordinate $distance" >> ../movable_atom.txt
            break 
        elif [ $k -eq 5 ]
        then
            sed -i "$line_number"'s/$/& F F F/' POSCAR
        fi
    done    
    done
    cd ..
done