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
    
    if np.abs(apos[2]) <= (radius):
      todelete.append(atom)


  xlength = (maxx-minx)/10
  zlength = (maxz-minz)/10.
  ylength = (maxy-miny)/10.

  porearea = 2*xlength*ylength
  numoh = np.floor(porearea*ohdensity)
  
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
    distcenter = apos[2]*apos[2]
    
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

  import pdb

  pdb.set_trace()

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
        angle = np.abs(apos[2])/apos[2]
        # add atom 1 \AA away from the oxygen, pointing to the central axis
        a = mol.NewAtom()
        a.SetAtomicNum(1) # hydrogen atom
        a.SetVector(apos[0], apos[1], apos[2] - 1*angle) # coordinates
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
