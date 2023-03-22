#!/bin/bash
V=64.56
step=17
rm polar_tot.dat
for((i=0;i<=$step;i++))
do
landa=`echo "scale=6;0.0+1.0/($step)*($i)"|bc` 
echo $landa >> ../landa.dat
wait
done

nstep=`expr $step + 1`
for((i=1;i<=$nstep;i++))
do
I1=$(sed -n "$i"p"" polar_ion.dat | awk '{printf  $1 }')
I2=$(sed -n "$i"p"" polar_ion.dat | awk '{printf  $2 }')
I3=$(sed -n "$i"p"" polar_ion.dat | awk '{printf  $3 }')
E1=$(sed -n "$i"p"" polar_ele.dat | awk '{printf  $1 }')
E2=$(sed -n "$i"p"" polar_ele.dat | awk '{printf  $2 }')
E3=$(sed -n "$i"p"" polar_ele.dat | awk '{printf  $3 }')
T1=`echo "scale=6;(($I1)-($E1))*1602.17/($V)"|bc` 
T2=`echo "scale=6;(($I2)-($E2))*1602.17/($V)"|bc` 
T3=`echo "scale=6;(($I3)-($E3))*1602.17/($V)"|bc` 
echo $T1 $T2 $T3  >> ./polar_tot.dat
wait
wait
done

#FILENAME="ploar_ion.dat"
#cat $FILENAME | while read LINE
#do
#I1=$(echo "$LINE" | awk '{printf  $1 }')
#I2=$(echo "$LINE" | awk '{printf  $2 }')
#I3=$(echo "$LINE" | awk '{printf  $3 }')
#echo $I1 $I2 $I3 >> ./test.dat
#done
#D=$(sed -n "3p" test.dat)
#echo $D
