#!/bin/bash
XYZ_FILE='sio2_cristobalite.xyz'

for j in 0 1 2 3 4;do
for i in 4 7.5 11 12.5 15 17.5 20;do 
    python cut_pore.py files/$XYZ_FILE $i -o pore_${i}_silanol_${j}.pdb --option plateh --silanol-density ${j} --buffer-size 2
    python cut_pore.py files/$XYZ_FILE $i -o pore_${i}_silanol_${j}.xyz --option plateh --silanol-density ${j} --buffer-size 2
    python rename_pdb_plate.py pore_${i}_silanol_${j}.pdb ${i} --detect-ionized-o --rename-label
    grep -v 'CONECT' rpore_${i}_silanol_${j}.pdb > pore_${i}_silanol_${j}.pdb

    mkdir -p outputs/pore_${i}_silanol_${j}/
    mv pore_${i}_silanol_${j}.xyz outputs/pore_${i}_silanol_${j}/
    mv pore_${i}_silanol_${j}.pdb outputs/pore_${i}_silanol_${j}/
    rm *.pdb

done
done 
