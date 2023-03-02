#!/bin/bash
XYZ_FILE='sio2_cristobalite.xyz'
#XYZ_FILE='sio2_amorphous.xyz'
buffer_size=3

#i=17.5
#j=4
for j in 0 1 2 3 4;do 
for i in 10 7.5 10 12.5 15 17.5 21;do 
#for i in 5 7.5 10 12.5 15 17.5 20;do 
    echo "Densidade silanol : $j" 
    echo "Distância da origem à parede do poro : $i"    
    python cut_pore.py files/$XYZ_FILE $i -o pore_${i}_silanol_${j}.pdb --option plateh --silanol-density ${j} --buffer-size $buffer_size
    python cut_pore.py files/$XYZ_FILE $i -o pore_${i}_silanol_${j}.xyz --option plateh --silanol-density ${j} --buffer-size $buffer_size
    python rename_pdb_plate.py pore_${i}_silanol_${j}.pdb ${i} --detect-ionized-o --rename-label --buffer-radius $buffer_size
    grep -v 'CONECT' rpore_${i}_silanol_${j}.pdb > pore_${i}_silanol_${j}.pdb
    mkdir -p outputs/pore_${i}_silanol_${j}/
    mv pore_${i}_silanol_${j}.xyz outputs/pore_${i}_silanol_${j}/
    mv pore_${i}_silanol_${j}.pdb outputs/pore_${i}_silanol_${j}/
    rm *.pdb
    echo ''
done
done 
