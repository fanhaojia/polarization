#!/bin/bash
V=64.31
rm polar_ion.dat polar_tot.dat polar_ele.dat
step=17
for((i=0;i<=17;i++))
do
cd $i
wait
I1=`grep 'Ionic dipole moment' OUTCAR | awk '{printf  $5 }'` 
I2=`grep 'Ionic dipole moment' OUTCAR | awk '{printf  $6 }'`
I3=`grep 'Ionic dipole moment' OUTCAR | awk '{printf  $7 }'`
echo $I1 $I2 $I3  >> ../polar_ion.dat
wait
E1=`grep 'Total electronic dipole moment' OUTCAR | awk '{printf  $6 }'` 
E2=`grep 'Total electronic dipole moment' OUTCAR | awk '{printf  $7 }'`
E3=`grep 'Total electronic dipole moment' OUTCAR | awk '{printf  $8 }'`
echo $E1 $E2 $E3  >> ../polar_ele.dat
wait
T1=`echo "scale=6;(($I1)-($E1))*1602.17/($V)"|bc` 
T2=`echo "scale=6;(($I2)-($E2))*1602.17/($V)"|bc` 
T3=`echo "scale=6;(($I3)-($E3))*1602.17/($V)"|bc` 
echo $T1 $T2 $T3  >> ../polar_tot.dat
wait
landa=`echo "scale=6;0.0+1.0/($step)*($i)"|bc` 
echo $landa >> ../landa.dat
wait
E=`grep "TOTEN" OUTCAR | tail -1 | awk  '{printf  "%12.6f \n" ,  $5 }'`
echo $E  >>../energy.dat
wait
#rm CHG* W* DOSCAR EIG* X* PC* RE* OS* CON* vasprun* PROCAR IB* ploar*
cd ../
done


