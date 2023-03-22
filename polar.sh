#!/bin/bash
step=17
for((i=0;i<=$step;i++))
do
cat >> KPOINTS <<!
888
0
Gamma
8 8 8
0 0 0
!
cat >> INCAR <<!
ISMEAR = 0; SIGMA = 0.05
LORBIT=11
ALGO = N;
ENCUT = 600
NELM = 170
EDIFF = 1E-6
NWRITE=3
NSW = 0 ;   POTIM=0.2 ; IBRION = 2 ; NFREE = 10 ;
################################################
LPEAD=.TRUE.
IPEAD= 4
LCALCPOL = .TRUE.
#LEFG=.TRUE.
##############################################
!
mkdir $i
cp ./PATH/CONTCAR$i.vasp ./$i/POSCAR
cp INCAR POTCAR KPONITS ./$i
cd $i
yhrun -p sz-renwei -N 3 -n 72 /THFS/home/renwei01/softwares/vasp/bin/vasp_std >report&
wait
rm CHG* W* DOSCAR EIG* X* PC* RE* OS* CON* vasprun* PROCAR IB* ploar*
cd ../
done
