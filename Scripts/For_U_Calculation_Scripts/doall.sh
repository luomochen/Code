#
# To run VASP this script calls $vasp_std
# (or posibly $vasp_gam and/or $vasp_ncl).
# These variables can be defined by sourcing vaspcmd
. vaspcmd 2> /dev/null

#
# When vaspcmd is not available and $vasp_std,
# $vasp_gam, and/or $vasp_ncl are not set as environment
# variables, you can specify them here
[ -z "`echo $vasp_std`" ] && vasp_std="mpirun -np 8 /opt/VTST_vasp.5.4.4/bin/vasp_std"
[ -z "`echo $vasp_gam`" ] && vasp_gam="mpirun -np 8 /opt/VTST_vasp.5.4.4/bin/vasp_gam"
[ -z "`echo $vasp_ncl`" ] && vasp_ncl="mpirun -np 8 /opt/VTST_vasp.5.4.4/bin/vasp_ncl"

#
# The real work starts here
#

cat > INCAR.DFT <<!
SYSTEM       = NiO AFM 

PREC         = A

EDIFF        = 1E-6
#AMIX         = 0.2
#BMIX         = 0.000001
#AMIX_MAG     = 0.2
#BMIX_MAG     = 0.000001
#NELM         = 150

ISMEAR       = 0
SIGMA        = 0.2

ISPIN        = 2
MAGMOM       = 2.0 -1.0  1.0 -1.0  1.0 \\
              -1.0  1.0 -1.0  1.0 -1.0 \\
               1.0 -1.0  1.0 -1.0  1.0 \\
              -1.0  1.0 -1.0  1.0 -1.0 \\
            16*0.0

LORBIT       = 11

LMAXMIX      = 4 
!

cp INCAR.DFT INCAR
rm WAVECAR CHGCAR

$vasp_std

cp OUTCAR  OUTCAR.0
cp OSZICAR OSZICAR.0
cp WAVECAR WAVECAR.0
cp CHGCAR  CHGCAR.0


for v in +0.05 -0.05 +0.10 -0.10 +0.15 -0.15 +0.20 -0.20
do

cp INCAR.DFT INCAR
cat >> INCAR <<!
ICHARG       = 11

LDAU         = .TRUE.
LDAUTYPE     =  3
LDAUL        =  2 -1 -1
LDAUU        =  $v 0.00 0.00
LDAUJ        =  $v 0.00 0.00
LDAUPRINT    =  2
!

cp WAVECAR.0 WAVECAR
cp CHGCAR.0  CHGCAR

$vasp_std

cp OSZICAR OSZICAR.V=$v.ICHARG=11
cp OUTCAR  OUTCAR.V=$v.ICHARG=11

cp INCAR.DFT INCAR
cat >> INCAR <<!
LDAU         = .TRUE.
LDAUTYPE     =  3
LDAUL        =  2 -1 -1
LDAUU        =  $v 0.00 0.00
LDAUJ        =  $v 0.00 0.00
LDAUPRINT    =  2
!

$vasp_std

cp OSZICAR OSZICAR.V=$v
cp OUTCAR  OUTCAR.V=$v

done
