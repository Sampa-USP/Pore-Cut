#!/usr/bin/env python3
"""
Receives a .xyz with amorphous silica and cut a
cylindrical pore with a given radius.

The pore is centered at 0.0 to be compatible with Cassandra.

Author: Henrique Musseli Cezar
Date: OCT/2021
"""

import sys
import argparse
import os
import numpy as np
from collections import OrderedDict
from openbabel import pybel
from openbabel import openbabel


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


def cut_pore(infile, radius, buffer, ohdensity, outfile):
  # set openbabel file format
  base, ext = os.path.splitext(infile)
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:], os.path.splitext(outfile)[1][1:])
  # trick to disable ring perception and make the ReadFile waaaay faster
  # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
  obConversion.AddOption("b", openbabel.OBConversion.INOPTIONS) 

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, infile)
  if ext[1:].lower() == "pdb":
    mol.SetFlag(openbabel.OB_PERIODIC_MOL) # makes the system periodic
  mol.ConnectTheDots() # necessary because of the 'b' INOPTION    
  mol.Center()

  # make a check of how many undercoordinate atoms there are in the input xyz
  if ext[1:].lower() == "pdb":
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Initial sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))

  # for each atom, check if it is in the cylinder
  todelete = []
  maxz = sys.float_info.min
  minz = sys.float_info.max
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    if apos[2] > maxz:
      maxz = apos[2]
    if apos[2] < minz:
      minz = apos[2]

    if (((apos[0]*apos[0]) + (apos[1]*apos[1])) <= (radius*radius)):
      todelete.append(atom)

  zlength = (maxz-minz)/10.
  porearea = 2.*np.pi*(radius/10.)*zlength
  numoh = np.floor(porearea*ohdensity)
  print("Estimated z length (nm): {}".format(zlength))
  print("Estimated pore area (nm^2): {}".format(porearea))
  print("Number of needed OH: {}".format(numoh))
  numoh = numoh - (numoh%4)
  print("Will add {} hydrogen atoms".format(numoh))
  print("Silanol density (nm^(-2)): {}".format(numoh/porearea))

  # delete atoms
  deleted = {8: 0, 14: 0}
  for atom in todelete:
    deleted[atom.GetAtomicNum()] += 1
    mol.DeleteAtom(atom)

  print("Deleted {} oxygen and {} silicon atoms.".format(deleted[8], deleted[14]))

  # add buffer to get the atoms at the surface
  rbuf = radius + buffer
  inbuffer = {}
  buffatoms = {8: 0, 14: 0}
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    distcenter = (apos[0]*apos[0]) + (apos[1]*apos[1])
    if (distcenter <= (rbuf*rbuf)):
      buffatoms[atom.GetAtomicNum()] += 1
      inbuffer[distcenter] = atom

  # create dictionary ordered by distance from pore center
  ordbuffer = OrderedDict(sorted(inbuffer.items()))
  # print(ordbuffer.keys())

  # delete first the O atoms with valence less than 1
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 1):
      deleted[8] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # delete first the Si atoms with valence less than 3
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 3):
      deleted[14] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[14] -= 1
    ordbuffer.pop(dist)

  # check if the number of oxygen atoms is enough and balance stechiometry
  # still need to delete some O atoms
  todelete = 2*deleted[14]-deleted[8]
  deleteddists = []
  if (todelete) > 0:
    excsilicon = 0
    ndeleted = 0
    # delete first the O atoms with valence less than 2
    for (dist, atom) in ordbuffer.items():
      if ndeleted >= todelete:
        break

      if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        deleted[8] += 1
        ndeleted += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

    # if we still need to delete more O, then delete the closest to the center
    if ndeleted < todelete:
      # delete atoms from buffer
      for dist in deleteddists:
        buffatoms[8] -= 1
        ordbuffer.pop(dist)
      deleteddists = []
      for (dist, atom) in ordbuffer.items():
        if ndeleted >= todelete:
          break

        if atom.GetAtomicNum() == 8:
          deleted[8] += 1
          ndeleted += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)

  # will delete more Si in the next step
  elif (todelete) < 0:
    # we must have an even number of O
    if todelete%2 == 1:
      deletedodd = False
      for (dist, atom) in ordbuffer.items():
        if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
          deletedodd = True
          deleted[8] += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)
          break

      if not deletedodd:
        for (dist, atom) in ordbuffer.items():
          if atom.GetAtomicNum() == 8:
            deletedodd = True
            deleted[8] += 1
            deleteddists.append(dist)
            mol.DeleteAtom(atom)
            break

    excsilicon = np.ceil(-todelete/2)
  # just delete the Si atoms in the next step
  else:
    excsilicon = 0

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # print(ordbuffer.keys())
  print("Deleted {} oxygen and {} silicon atoms after fix.".format(deleted[8], deleted[14]))
  print("Excess silicon atoms: {}".format(excsilicon))

  # for each Si we remove, we add 2 silanol
  # first remove the Si atoms coordinated with less than 4 atoms
  numsilicondelete = np.floor(numoh/4) + excsilicon
  print("Will delete {} silicon atoms".format(numsilicondelete))
  print("And add {} hydrogen atoms".format(numoh))
  ndelsi = 0
  deleteddists = []
  for (dist, atom) in ordbuffer.items():
    if ndelsi >= numsilicondelete:
      break

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
      ndelsi += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # if we still need to delete Si atoms, delete them
  deleteddists = []
  if ndelsi < numsilicondelete:
    for (dist, atom) in ordbuffer.items():
      if ndelsi >= numsilicondelete:
        break

      if atom.GetAtomicNum() == 14:
        ndelsi += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # now add the H connected to the more central oxygens
  # also checks if the O is not connected to an Si that already has an O-H
  naddedh = 0
  silist = []
  for (dist, atom) in ordbuffer.items():
    if naddedh >= numoh:
      break

    # get the bonded atom (if it's an O it should only be bonded to 1 atom to become an O-H)
    if ((atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() == 1)):
      for atom2 in openbabel.OBAtomAtomIter(atom):
        nbrAtom = atom2

      if nbrAtom.GetId() not in silist:
        # add the bonded Si atom to list
        if nbrAtom.GetAtomicNum() == 14:
          silist.append(nbrAtom.GetId())

        # get the angle of atom (cylindrical coordinates)
        apos = (atom.GetX(), atom.GetY(), atom.GetZ())
        angle = np.arctan2(apos[1], apos[0])
        # add atom 1 \AA away from the oxygen, pointing to the central axis
        a = mol.NewAtom()
        a.SetAtomicNum(1) # hydrogen atom
        a.SetVector(apos[0]-np.cos(angle), apos[1]-np.sin(angle), apos[2]) # coordinates
        naddedh += 1

    else:
      continue

  print("Number of added hydrogens: {}".format(naddedh))

  # make a check of how many undercoordinate atoms there are in the output xyz
  if ext[1:].lower() == "pdb":
    mol.ConnectTheDots() # necessary because of deleted atoms    
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Final sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))


  # write final structure
  obConversion.WriteFile(mol, outfile)
  
  
def cut_plate(infile, radius, buffer, ohdensity, outfile):
  # set openbabel file format
  base, ext = os.path.splitext(infile)
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:], os.path.splitext(outfile)[1][1:])
  # trick to disable ring perception and make the ReadFile waaaay faster
  # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
  obConversion.AddOption("b", openbabel.OBConversion.INOPTIONS) 

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, infile)
  if ext[1:].lower() == "pdb":
    mol.SetFlag(openbabel.OB_PERIODIC_MOL) # makes the system periodic
  mol.ConnectTheDots() # necessary because of the 'b' INOPTION    
  mol.Center()

  # make a check of how many undercoordinate atoms there are in the input xyz
  if ext[1:].lower() == "pdb":
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Initial sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))

  # for each atom, check if it is in the cylinder
  todelete = []
  maxz = sys.float_info.min
  minz = sys.float_info.max
  
  maxx = sys.float_info.min
  minx = sys.float_info.max
  
  maxy = sys.float_info.min
  miny = sys.float_info.max
  
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    if apos[2] > maxz:
      maxz = apos[2]
    if apos[2] < minz:
      minz = apos[2]

    if apos[0] > maxx:
      maxx = apos[0]
    if apos[0] < minx:
      minx = apos[0]

    if apos[1] > maxy:
      maxy = apos[1]
    if apos[1] < miny:
      miny = apos[1]
    
    if (((apos[2]*apos[2])) <= (radius*radius)):
      todelete.append(atom)


  xlength = (maxx-minx)/10
  zlength = (maxz-minz)/10.
  ylength = (maxy-miny)/10.
  porearea = 2*xlength*zlength
  numoh = np.floor(porearea*ohdensity)
  print("Estimated z length (nm): {}".format(zlength))
  print("Estimated y length (nm): {}".format(ylength))
  print("Estimated x length (nm): {}".format(xlength))
  
  print("Min X (nm): {}".format(minx))
  print("Max X (nm): {}".format(maxx))
  
  print("Mix Y (nm): {}".format(miny))
  print("Max Y (nm): {}".format(maxy))
  
  print("Min Z (nm): {}".format(minz))
  print("Max Z (nm): {}".format(maxz))
  
  print("Estimated pore area (nm^2): {}".format(porearea))
  print("Number of needed OH: {}".format(numoh))
  numoh = numoh - (numoh%4)
  print("Will add {} hydrogen atoms".format(numoh))
  print("Silanol density (nm^(-2)): {}".format(numoh/porearea))

  # delete atoms
  deleted = {8: 0, 14: 0}
  for atom in todelete:
    deleted[atom.GetAtomicNum()] += 1
    mol.DeleteAtom(atom)

  print("Deleted {} oxygen and {} silicon atoms.".format(deleted[8], deleted[14]))

  # add buffer to get the atoms at the surface
  rbuf = radius + buffer
  inbuffer = {}
  buffatoms = {8: 0, 14: 0}
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    distcenter = apos[1]*apos[1]
    
    if (distcenter <= (rbuf*rbuf)):
      buffatoms[atom.GetAtomicNum()] += 1
      inbuffer[distcenter] = atom

  # create dictionary ordered by distance from pore center
  ordbuffer = OrderedDict(sorted(inbuffer.items()))
  # print(ordbuffer.keys())

  # delete first the O atoms with valence less than 1
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 1):
      deleted[8] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # delete first the Si atoms with valence less than 3
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 3):
      deleted[14] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[14] -= 1
    ordbuffer.pop(dist)

  # check if the number of oxygen atoms is enough and balance stechiometry
  # still need to delete some O atoms
  todelete = 2*deleted[14]-deleted[8]
  deleteddists = []
  if (todelete) > 0:
    excsilicon = 0
    ndeleted = 0
    # delete first the O atoms with valence less than 2
    for (dist, atom) in ordbuffer.items():
      if ndeleted >= todelete:
        break

      if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        deleted[8] += 1
        ndeleted += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

    # if we still need to delete more O, then delete the closest to the center
    if ndeleted < todelete:
      # delete atoms from buffer
      for dist in deleteddists:
        buffatoms[8] -= 1
        ordbuffer.pop(dist)
      deleteddists = []
      for (dist, atom) in ordbuffer.items():
        if ndeleted >= todelete:
          break

        if atom.GetAtomicNum() == 8:
          deleted[8] += 1
          ndeleted += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)

  # will delete more Si in the next step
  elif (todelete) < 0:
    # we must have an even number of O
    if todelete%2 == 1:
      deletedodd = False
      for (dist, atom) in ordbuffer.items():
        if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
          deletedodd = True
          deleted[8] += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)
          break

      if not deletedodd:
        for (dist, atom) in ordbuffer.items():
          if atom.GetAtomicNum() == 8:
            deletedodd = True
            deleted[8] += 1
            deleteddists.append(dist)
            mol.DeleteAtom(atom)
            break

    excsilicon = np.ceil(-todelete/2)
  # just delete the Si atoms in the next step
  else:
    excsilicon = 0

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # print(ordbuffer.keys())
  print("Deleted {} oxygen and {} silicon atoms after fix.".format(deleted[8], deleted[14]))
  print("Excess silicon atoms: {}".format(excsilicon))

  # for each Si we remove, we add 2 silanol
  # first remove the Si atoms coordinated with less than 4 atoms
  numsilicondelete = np.floor(numoh/4) + excsilicon
  print("Will delete {} silicon atoms".format(numsilicondelete))
  print("And add {} hydrogen atoms".format(numoh))
  ndelsi = 0
  deleteddists = []
  for (dist, atom) in ordbuffer.items():
    if ndelsi >= numsilicondelete:
      break

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
      ndelsi += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # if we still need to delete Si atoms, delete them
  deleteddists = []
  if ndelsi < numsilicondelete:
    for (dist, atom) in ordbuffer.items():
      if ndelsi >= numsilicondelete:
        break

      if atom.GetAtomicNum() == 14:
        ndelsi += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # now add the H connected to the more central oxygens
  # also checks if the O is not connected to an Si that already has an O-H
  naddedh = 0
  silist = []
  for (dist, atom) in ordbuffer.items():
    if naddedh >= numoh:
      break

    # get the bonded atom (if it's an O it should only be bonded to 1 atom to become an O-H)
    if ((atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() == 1)):
      for atom2 in openbabel.OBAtomAtomIter(atom):
        nbrAtom = atom2

      if nbrAtom.GetId() not in silist:
        # add the bonded Si atom to list
        if nbrAtom.GetAtomicNum() == 14:
          silist.append(nbrAtom.GetId())

        # get the angle of atom (cylindrical coordinates)
        apos = (atom.GetX(), atom.GetY(), atom.GetZ())
        angle = np.arctan2(apos[1], apos[0])
        # add atom 1 \AA away from the oxygen, pointing to the central axis
        a = mol.NewAtom()
        a.SetAtomicNum(1) # hydrogen atom
        a.SetVector(apos[0]-np.cos(angle), apos[1]-np.sin(angle), apos[2]) # coordinates
        naddedh += 1

    else:
      continue

  print("Number of added hydrogens: {}".format(naddedh))

  # make a check of how many undercoordinate atoms there are in the output xyz
  if ext[1:].lower() == "pdb":
    mol.ConnectTheDots() # necessary because of deleted atoms    
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Final sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))


  # write final structure
  obConversion.WriteFile(mol, outfile)


def cut_plate_redirect_o2(infile, radius, buffer, ohdensity, outfile):
  # set openbabel file format
  base, ext = os.path.splitext(infile)
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:], os.path.splitext(outfile)[1][1:])
  # trick to disable ring perception and make the ReadFile waaaay faster
  # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
  obConversion.AddOption("b", openbabel.OBConversion.INOPTIONS) 

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, infile)
  if ext[1:].lower() == "pdb":
    mol.SetFlag(openbabel.OB_PERIODIC_MOL) # makes the system periodic
  mol.ConnectTheDots() # necessary because of the 'b' INOPTION    
  mol.Center()

  # make a check of how many undercoordinate atoms there are in the input xyz
  if ext[1:].lower() == "pdb":
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Initial sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))

  # for each atom, check if it is in the cylinder
  todelete = []
  maxz = sys.float_info.min
  minz = sys.float_info.max
  
  maxx = sys.float_info.min
  minx = sys.float_info.max
  
  maxy = sys.float_info.min
  miny = sys.float_info.max
  
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    if apos[2] > maxz:
      maxz = apos[2]
    if apos[2] < minz:
      minz = apos[2]

    if apos[0] > maxx:
      maxx = apos[0]
    if apos[0] < minx:
      minx = apos[0]

    if apos[1] > maxy:
      maxy = apos[1]
    if apos[1] < miny:
      miny = apos[1]
    
    if (((apos[2]*apos[2])) <= (radius*radius)):
      todelete.append(atom)


  xlength = (maxx-minx)/10
  zlength = (maxz-minz)/10.
  ylength = (maxy-miny)/10.
  porearea = 2*xlength*zlength
  numoh = np.floor(porearea*ohdensity)
  print("Estimated z length (nm): {}".format(zlength))
  print("Estimated y length (nm): {}".format(ylength))
  print("Estimated x length (nm): {}".format(xlength))
  
  print("Min X (nm): {}".format(minx))
  print("Max X (nm): {}".format(maxx))
  
  print("Mix Y (nm): {}".format(miny))
  print("Max Y (nm): {}".format(maxy))
  
  print("Min Z (nm): {}".format(minz))
  print("Max Z (nm): {}".format(maxz))
  
  print("Estimated pore area (nm^2): {}".format(porearea))
  print("Number of needed OH: {}".format(numoh))
  numoh = numoh - (numoh%4)
  print("Will add {} hydrogen atoms".format(numoh))
  print("Silanol density (nm^(-2)): {}".format(numoh/porearea))

  # delete atoms
  deleted = {8: 0, 14: 0}
  for atom in todelete:
    deleted[atom.GetAtomicNum()] += 1
    mol.DeleteAtom(atom)

  print("Deleted {} oxygen and {} silicon atoms.".format(deleted[8], deleted[14]))

  # add buffer to get the atoms at the surface
  rbuf = radius + buffer
  inbuffer = {}
  buffatoms = {8: 0, 14: 0}
  for atom in openbabel.OBMolAtomIter(mol):
    apos = (atom.GetX(), atom.GetY(), atom.GetZ())
    distcenter = apos[1]*apos[1]
    
    if (distcenter <= (rbuf*rbuf)):
      buffatoms[atom.GetAtomicNum()] += 1
      inbuffer[distcenter] = atom

  # create dictionary ordered by distance from pore center
  ordbuffer = OrderedDict(sorted(inbuffer.items()))
  # print(ordbuffer.keys())

  # delete first the O atoms with valence less than 1
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 1):
      deleted[8] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # delete first the Si atoms with valence less than 3
  deleteddists = []
  for (dist, atom) in ordbuffer.items():

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 3):
      deleted[14] += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[14] -= 1
    ordbuffer.pop(dist)

  # check if the number of oxygen atoms is enough and balance stechiometry
  # still need to delete some O atoms
  todelete = 2*deleted[14]-deleted[8]
  deleteddists = []
  if (todelete) > 0:
    excsilicon = 0
    ndeleted = 0
    # delete first the O atoms with valence less than 2
    for (dist, atom) in ordbuffer.items():
      if ndeleted >= todelete:
        break

      if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        deleted[8] += 1
        ndeleted += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

    # if we still need to delete more O, then delete the closest to the center
    if ndeleted < todelete:
      # delete atoms from buffer
      for dist in deleteddists:
        buffatoms[8] -= 1
        ordbuffer.pop(dist)
      deleteddists = []
      for (dist, atom) in ordbuffer.items():
        if ndeleted >= todelete:
          break

        if atom.GetAtomicNum() == 8:
          deleted[8] += 1
          ndeleted += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)

  # will delete more Si in the next step
  elif (todelete) < 0:
    # we must have an even number of O
    if todelete%2 == 1:
      deletedodd = False
      for (dist, atom) in ordbuffer.items():
        if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
          deletedodd = True
          deleted[8] += 1
          deleteddists.append(dist)
          mol.DeleteAtom(atom)
          break

      if not deletedodd:
        for (dist, atom) in ordbuffer.items():
          if atom.GetAtomicNum() == 8:
            deletedodd = True
            deleted[8] += 1
            deleteddists.append(dist)
            mol.DeleteAtom(atom)
            break

    excsilicon = np.ceil(-todelete/2)
  # just delete the Si atoms in the next step
  else:
    excsilicon = 0

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # print(ordbuffer.keys())
  print("Deleted {} oxygen and {} silicon atoms after fix.".format(deleted[8], deleted[14]))
  print("Excess silicon atoms: {}".format(excsilicon))

  # for each Si we remove, we add 2 silanol
  # first remove the Si atoms coordinated with less than 4 atoms
  numsilicondelete = np.floor(numoh/4) + excsilicon
  print("Will delete {} silicon atoms".format(numsilicondelete))
  print("And add {} hydrogen atoms".format(numoh))
  ndelsi = 0
  deleteddists = []
  for (dist, atom) in ordbuffer.items():
    if ndelsi >= numsilicondelete:
      break

    if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
      ndelsi += 1
      deleteddists.append(dist)
      mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # if we still need to delete Si atoms, delete them
  deleteddists = []
  if ndelsi < numsilicondelete:
    for (dist, atom) in ordbuffer.items():
      if ndelsi >= numsilicondelete:
        break

      if atom.GetAtomicNum() == 14:
        ndelsi += 1
        deleteddists.append(dist)
        mol.DeleteAtom(atom)

  # delete atoms from buffer
  for dist in deleteddists:
    buffatoms[8] -= 1
    ordbuffer.pop(dist)

  # now add the H connected to the more central oxygens
  # also checks if the O is not connected to an Si that already has an O-H
  naddedh = 0
  silist = []
  for (dist, atom) in ordbuffer.items():
    if naddedh >= numoh:
      #break
      pass
    # get the bonded atom (if it's an O it should only be bonded to 1 atom to become an O-H)
    if ((atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() == 1)):
      for atom2 in openbabel.OBAtomAtomIter(atom):
        nbrAtom = atom2

      if nbrAtom.GetId() not in silist:
        # add the bonded Si atom to list
        if nbrAtom.GetAtomicNum() == 14:
          silist.append(nbrAtom.GetId())

        # get the angle of atom (cylindrical coordinates)
        apos = (atom.GetX(), atom.GetY(), atom.GetZ())
        angle = np.arctan2(apos[1], apos[0])
        # add atom 1 \AA away from the oxygen, pointing to the central axis
        a = mol.NewAtom()
        a.SetAtomicNum(8) # hydrogen atom
        a.SetVector(apos[0]-np.cos(angle), apos[1]-np.sin(angle), apos[2]) # coordinates
        naddedh += 1

    else:
      continue

  print("Number of added hydrogens: {}".format(naddedh))

  # make a check of how many undercoordinate atoms there are in the output xyz
  if ext[1:].lower() == "pdb":
    mol.ConnectTheDots() # necessary because of deleted atoms    
    uSi = 0
    uO = 0
    for atom in openbabel.OBMolAtomIter(mol):
      if (atom.GetAtomicNum() == 14) and (atom.GetExplicitDegree() < 4):
        uSi += 1
      elif (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() < 2):
        uO += 1
    print("Final sample has {} under-coordinated Si and {} under-coordinated O.".format(uSi, uO))


  # write final structure
  obConversion.WriteFile(mol, outfile)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Receives a .xyz or file read by OpenBabel and cut a cylindrical pore with a given radius at 0.0")
  parser.add_argument("xyzfile", type=extant_file, help="path to a .xyz containing amorphous silica")
  parser.add_argument("radius", type=float, help="radius of the pore")
  parser.add_argument("-o", "--output", metavar="OUTFILE" , help="name of output file (default = pore.pdb)", default="pore.pdb")
  parser.add_argument("--buffer-size", type=float, help="buffer added to radius to get the surface atoms (default = 1.0)", default=1.0)
  parser.add_argument("--silanol-density", type=float, help="density of Si-O-H terminations at the pore surface in nm^(-2) (default = 3)", default=100)
  parser.add_argument("--option", type=str, help="options are : cylinder_oh_ | plate_oh_ | plate_o_",default="plateo")

  args = parser.parse_args()

  obabel_sup = ["acr", "adf", "adfout", "alc", "arc", "bgf", "box", "bs", "c3d1", "c3d2", "cac", "caccrt", "cache", "cacint", "can", "car", "ccc", "cdx", "cdxml", "cht", "cif", "ck", "cml", "cmlr", "com", "copy", "crk2d", "crk3d", "csr", "cssr", "ct", "cub", "cube", "dmol", "dx", "ent", "fa", "fasta", "fch", "fchk", "fck", "feat", "fh", "fix", "fpt", "fract", "fs", "fsa", "g03", "g92", "g94", "g98", "gal", "gam", "gamin", "gamout", "gau", "gjc", "gjf", "gpr", "gr96", "gukin", "gukout", "gzmat", "hin", "inchi", "inp", "ins", "jin", "jout", "mcdl", "mcif", "mdl", "ml2", "mmcif", "mmd", "mmod", "mol", "mol2", "molden", "molreport", "moo", "mop", "mopcrt", "mopin", "mopout", "mpc", "mpd", "mpqc", "mpqcin", "msi", "msms", "nw", "nwo", "outmol", "pc", "pcm", "pdb", "png", "pov", "pqr", "pqs", "prep", "qcin", "qcout", "report", "res", "rsmi", "rxn", "sd", "sdf", "smi", "smiles", "sy2", "t41", "tdd", "test", "therm", "tmol", "txt", "txyz", "unixyz", "vmol", "xed", "xml", "xyz", "yob", "zin"]

  # get basename and file extension
  base, ext = os.path.splitext(args.xyzfile)

  print(args)


  if ext[1:] not in obabel_sup:
    sys.exit("Error: {} files are not accepted by OpenBabel.".format(ext[1:]))

  if args.option=="poreh":
    cut_pore(args.xyzfile, args.radius, args.buffer_size, args.silanol_density, str(args.radius) +"_" + str(args.silanol_density) + "_" +args.output)
  elif args.option=="plateh":
    cut_plate(args.xyzfile, args.radius, args.buffer_size, args.silanol_density,str(args.radius) +"_" + str(args.silanol_density) + "_" +args.output)
  elif args.option=="plateo":
    cut_plate_redirect_o2(args.xyzfile, args.radius, args.buffer_size, args.silanol_density,str(args.radius) +"_" + str(args.silanol_density) + "_" +args.output)

  else:
    sys.exit("Error! Neither @pore@ or @plate@ options has been selected.You need to choose one of them.")
