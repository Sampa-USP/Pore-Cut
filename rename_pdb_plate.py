"""
Receives a pdb and renames the Oxygens that are bonded to hydrogens.

Author: Henrique Musseli Cezar
Date: OCT/2021
"""

import os
import argparse
import sys
import numpy as np
try:
    import pybel
    import openbabel
except:
    from openbabel import pybel
    from openbabel import openbabel
import timeit

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


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a pdb and renames the oxygens that are bonded to hydrogens")
  parser.add_argument("pdbfile", type=extant_file, help="pdb containing the structure")
  parser.add_argument("radius", type=float, help="pore radius in Angstroms to consider the surface")
  parser.add_argument("--detect-ionized-o", action="store_true", help="detect ionized oxygens on the surface of the pore")  
  parser.add_argument("--buffer-radius", type=float, help="buffer added to pore radius to find surface atoms (default = 1.0)", default=1.0)
  parser.add_argument("--rename-label", action="store_true", help="rename the label instead of atom name")  


  args = parser.parse_args()

  rbufsq = (args.radius+args.buffer_radius)*(args.radius+args.buffer_radius)

   # set openbabel file format
  base, ext = os.path.splitext(args.pdbfile)
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:],"pdb")
  # trick to disable ring perception and make the ReadFile waaaay faster
  # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
  obConversion.AddOption("b", openbabel.OBConversion.INOPTIONS) 
  obConversion.AddOption("c", openbabel.OBConversion.INOPTIONS) # ignore the CONECTs

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, args.pdbfile)
  #mol.SetFlag(openbabel.OB_PERIODIC_MOL)


  mol.ConnectTheDots() # necessary because of the 'b' INOPTION

  # for each bond, check if it is an hydrogen connected to an oxygen
  oxygensh = []
  for bond in openbabel.OBMolBondIter(mol):
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    b1 = a1.GetAtomicNum()
    b2 = a2.GetAtomicNum()

    if ((b1 == 1) and (b2 == 8)):
      oxygensh.append(a2.GetId()+1)
    elif (b1 == 8) and (b2 == 1):
      oxygensh.append(a1.GetId()+1)

  # for each oxygen, check if it's an ionized oxygen of a silanol at the pore surface
  ionoxygen = []
  ionsilicon = []
  if args.detect_ionized_o:
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() == 1):
        z = atom.GetZ()
        if ((z*z) <= rbufsq):
          ionoxygen.append(atom.GetId()+1)
          for neigh in openbabel.OBAtomAtomIter(atom):
            ionsilicon.append(neigh.GetId()+1)

  print("Found {} OH, {} OI and {} SiI atoms".format(len(oxygensh), len(ionoxygen), len(ionsilicon)))

  # now rename based on the list
  anum = 0
  fout = open('r'+args.pdbfile, "w")
  with open(args.pdbfile, "r") as f:
    for line in f:
      if "HETATM" in line or "ATOM" in line:
        anum += 1
        if anum in oxygensh:
          llist = list(line.rstrip())
          if args.rename_label:
            llist[12] = "O"
            llist[13] = "H"
            line = "".join(llist)
          else:
            line = "".join(llist)+" OH"
        elif anum in ionoxygen:
          llist = list(line.rstrip())
          if args.rename_label:
            llist[12] = "O"
            llist[13] = "I"
            line = "".join(llist)
          else:
            line = "".join(llist)+" OI"
        elif anum in ionsilicon:
          llist = list(line.rstrip())
          if args.rename_label:
            llist[12] = "S"
            llist[13] = "i"
            llist[14] = "I"
            line = "".join(llist)
          else:
            line = "".join(llist)+" SiI"
        else:
          if args.rename_label:
            line = line.rstrip()
          else:
            line = line.rstrip() + " %s" % (line.split()[-1])
        fout.write(line+"\n")
      elif "CONECT" in line:
        pass
      else:
        fout.write(line)
