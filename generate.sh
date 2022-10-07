#!/bin/bash

for i in 5 7.5 10 12.5 15 17.5 20;do 
    python cut_pore.py files/sio2.xyz $i -o pore_${i}_silanol_5.pdb --option plateh --silanol-density 5
    python cut_pore.py files/sio2.xyz $i -o pore_${i}_silanol_5.xyz --option plateh --silanol-density 5
    python rename_pdb_plate.py pore_${i}_silanol_5.pdb ${i} --detect-ionized-o --rename-label
done

mv *.xyz output/
mv *.pdb output/
