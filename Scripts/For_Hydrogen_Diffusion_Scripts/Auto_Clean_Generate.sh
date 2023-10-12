#!/bin/bash
#----------------------------------------------------------
# 该脚本用于自动清理neb计算中产生的部分无用文件。
# 同时将重新生成一个用于计算的文件夹。
#----------------------------------------------------------

# 备份上一次计算。
cp -r path2 path2.bak
# 清理文件夹中的无用文件。
cd path2||exit
rm -rf exts.dat mep.eps movie movie.vasp movie.xyz neb.dat nebef.dat spline.dat vaspgr vasprun.xml vasp_test.*
list="01 02 03 04 05"
for i in $list
do
    cd "$i"||exit
    rm -rf CHG CHGCAR DOSCAR EIGENVAL fe.dat IBZKPT OSZICAR OUTCAR PCDAT POSCAR POSCAR.vasp POSCAR.xyz REPORT stdout WAVECAR XDATCAR
    mv CONTCAR POSCAR
    cd ..
done