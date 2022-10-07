"""
Receives a pdb and renames the Oxygens that are bonded to hydrogens.

Author: Henrique Musseli Cezar
Date: MAR/2020
"""

import os
import argparse
import sys
import numpy as np
try:
  import pybel
except:
  from openbabel import pybel

# from https://stackoverflow.com/a/11541495
def extant_file(x):
  """
  'Type' for argparse - checks that file exists but does not open.
  """
  if not os.path.exists(x):
      # Argparse uses the ArgumentTypeError to give a rejection message like:
      # error: argument input: x does not exist
      raise argparse.ArgumentTypeError("{0} does not exist".format(x))
  return x

def readxyz(file):
  species = []
  atoms = []
  with open(file,"r") as f:
    try:
      natoms = int(f.readline())
    except:
      print("The first line of the xyz file should be the number of atoms")
      sys.exit(0)

    # skip comment and read atomic positions and species
    f.readline()
    for line in f:
      species.append(line.split()[0])
      atoms.append(np.array([float(x) for x in line.split()[1:]]))

  return species, np.array(atoms)

def writexyz(species, atoms, file):
  with open(file,"w") as f:
    f.write("%d\n\n" % (len(species)))
    for sp, acoord in zip(species, atoms):
      f.write("%s\t%.6f\t%.6f\t%.6f\n" % (sp, acoord[0], acoord[1], acoord[2]))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a pdb and renames the oxygens that are bonded to hydrogens")
  parser.add_argument("pdbfile", type=extant_file, help="pdb containing the structure")
  parser.add_argument("--rename-label", action="store_true", help="rename the label instead of atom name")  

  args = parser.parse_args()


  import pdb
  pdb.set_trace()
  # read the pdb
  mol = pybel.readfile("pdb", args.pdbfile).__next__()

  ### HBD ### # Lipinski RO5 definition
  HBD_RO5 = pybel.Smarts("[#8;H1;!h1]") # Any oxygen with one attached hydrogen that is not implicit

  # get the oxygens
  oxygensh = HBD_RO5.findall(mol)
  oxygensh = [x[0] for x in oxygensh]

  # now rename based on the list
  anum = 0
  fout = open("labels_renamed.pdb", "w")
  with open(args.pdbfile, "r") as f:
    for line in f:
      if "HETATM" in line or "ATOM" in line:
        anum += 1
        if anum in oxygensh:
          llist = list(line)
          if args.rename_label:
            llist[12] = "O"
            llist[13] = "H"
          else:
            llist[76] = "O"
            llist[77] = "H"
          line = "".join(llist)
        fout.write(line)
      elif "CONECT" in line:
        pass
      else:
        fout.write(line)