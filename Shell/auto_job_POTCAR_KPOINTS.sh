#!/bin/bash

if [ "$potential_cut" == "n" ];then
   for ((i=0;i<"${#new_elements_type[*]}";i++));do
       potcar_elements="${new_elements_type[$i]}"
       #echo $potcar_elements
       if [ "$ex_co_po" == "GGA" ];then
           if [ $i == 0 ];then
              cat /home/xwm/VASP54_PAW_PBE/$potcar_elements/POTCAR > POTCAR
           else
              cat /home/xwm/VASP54_PAW_PBE/$potcar_elements/POTCAR >> POTCAR
           fi
        else
           if [ $i == 0 ];then
              cat /home/xwm/VASP_PAW_LDA/$potcar_elements/POTCAR > POTCAR
           else
              cat /home/xwm/VASP_PAW_LDA/$potcar_elements/POTCAR >> POTCAR
           fi
       fi
   done
elif [ "$potential_cut" == "y" ];then
   potcar_types=(`sed -n '3p' pre_POTCAR_KPOINTS`)
   #echo "${potcar_types[*]}"
   #echo "${#potcar_types[*]}"
   for ((i=0;i<"${#potcar_types[*]}";i++));do
       potcar_elements="${potcar_types[$i]}"
       if [ "$ex_co_po" == "GGA" ];then
           if [ $i == 0 ];then
              cat /home/xwm/VASP54_PAW_PBE/$potcar_elements/POTCAR > POTCAR
           else
              cat /home/xwm/VASP54_PAW_PBE/$potcar_elements/POTCAR >> POTCAR
           fi
        else
           if [ $i == 0 ];then
              cat /home/xwm/VASP_PAW_LDA/$potcar_elements/POTCAR > POTCAR
           else
              cat /home/xwm/VASP_PAW_LDA/$potcar_elements/POTCAR >> POTCAR
           fi
        fi

   done
else
   echo "You need to set 'potential_cut' and 'ex_co_po' to choose the exchange correlation potential" 
fi
mv POTCAR cal_U/relax


X=(`sed -n '3p' "$input_poscar"`)
Y=(`sed -n '4p' "$input_poscar"`)
Z=(`sed -n '5p' "$input_poscar"`)

x=(`echo "sqrt(${X[0]}*${X[0]}+${X[1]}*${X[1]}+${X[2]}*${X[2]})"|bc |awk '{printf "%.10f", $0}'`)
y=(`echo "sqrt(${Y[0]}*${Y[0]}+${Y[1]}*${Y[1]}+${Y[2]}*${Y[2]})"|bc |awk '{printf "%.10f", $0}'`)
z=(`echo "sqrt(${Z[0]}*${Z[0]}+${Z[1]}*${Z[1]}+${Z[2]}*${Z[2]})"|bc |awk '{printf "%.10f", $0}'`)

k_x=(`echo "40/$x"|bc | awk '{printf "%.0f", $0}'`)
k_y=(`echo "40/$y"|bc | awk '{printf "%.0f", $0}'`)
k_z=(`echo "40/$z"|bc | awk '{printf "%.0f", $0}'`)
kpoints=(`echo $k_x $k_y $k_z`)
sed -n '9,11'p pre_POTCAR_KPOINTS > KPOINTS
echo "${kpoints[*]}" >> KPOINTS
sed -n '13p' pre_POTCAR_KPOINTS >> KPOINTS
mv KPOINTS cal_U/relax 
mv IBZKPT cal_U/relax/KPOINTS              
