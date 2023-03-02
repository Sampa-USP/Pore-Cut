#!/bin/bash

main_file=$1

for file in `ls $main_file`;do 
for pdb in `ls $main_file/$file/*.pdb`;do

yes LJ | mcfgen.py $pdb --ffTemplate

ff_file="${pdb//pdb/ff}"
IFS='/' read -ra my_array<<<$ff_file
ff_file=${my_array[-1]}

IFS='_' read -ra my_array<<<$ff_file
size=${my_array[1]}

xyz_file="${pdb//pdb/xyz}"
IFS='/' read -ra my_array<<<$xyz_file
xyz_file=${my_array[-1]}

mcf_file="${pdb//pdb/mcf}"
IFS='/' read -ra my_array<<<$mcf_file
mcf_file=${my_array[-1]}


sed -z -i 's/Si\nSigma \nEpsilon \natom_type_charge /Si\nSigma 4.15\nEpsilon 0.389112\nnatom_type_charge +1.1/' $ff_file
sed -z -i 's/O\nSigma \nEpsilon \natom_type_charge /O\nSigma 3.47\nEpsilon 0.225936\nnatom_type_charge -0.55/' $ff_file

sed -z -i 's/SiI\nSigma \nEpsilon \natom_type_charge /SiI\nSigma 4.15\nEpsilon 0.389112\nnatom_type_charge +1.1/' $ff_file
sed -z -i 's/OI\nSigma \nEpsilon \natom_type_charge /OI\nSigma 3.47\nEpsilon 0.510448\nnatom_type_charge -0.675/' $ff_file
sed -z -i 's/H\nSigma \nEpsilon \natom_type_charge /H\nSigma 1.085\nEpsilon 0.06276\nnatom_type_charge +0.4/' $ff_file


echo $size
echo $xyz_file

sed "s/restricted_insertion slitpore 11/restricted_insertion slitpore $size/" gcmc2.inp > gcmc.inp
sed -z -i "s/gcmc_run_20A/gcmc_run_$size/" gcmc.inp
sed -z -i "s/plate_o_pore.xyz/$xyz_file/" gcmc.inp
sed -z -i "s/plate_o_pore.mcf/$mcf_file/" gcmc.inp


mcfgen.py $pdb
mv *.ff* $main_file/$file/
mv *.mcf $main_file/$file/
#rm $main_file/$file/gcmc2.inp

cp 2run.sh $main_file/$file/
cp gcmc.inp $main_file/$file/
cp spce.dat $main_file/$file/
cp waterSPCE $main_file/$file/waterSPCE.mcf

done
done
