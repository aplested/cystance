#!/usr/bin/env python

import __main__
import os
import pymol
import random
from pymol import cmd
import math
from datetime import datetime
from cysTools import formatLine, getAccess, getClashes, rearrange_subunits, hitCalculation

# Calling PyMOL in this way is done in the Anaconda Python distro
# bundled with Pymol 2.0 from Schroedinger
# /Applications/PyMOL.app/Contents/bin/python

# This script doesn't generate any new PDB models.
# It only measures the same properties as buryMyCys from existing PDB files
# PDB files for use with this script must be prepared!
# No segIDs! Sorry
# Chains A, B, C & D
##
# to get access to CCP4 the script should be run from the terminal with CCP4 sourced.
# use NCONT to count close contacts for this residue ("buried")
# use AREAIMOL to see if residue of interest (e.g. Cys 666) has been buried
# [proposed] use CNS to do simulated annealing on models and tidy them up (like in Morph)
# <minimize powell> is the command (from morph_dist.inp)
##
# Tell PyMOL to launch quiet (-q), and with no GUI (-c)
# Make sure PyMOL has finished to launch before using any PyMOL modules.
__main__.pymol_argv = [ 'pymol', '-qc' ]
pymol.finish_launching()

wd = os.path.expanduser('~') #working directory == home directory. Messy.
cmd.cd(wd)

#prepare logfile for results only
logfile = "log_" + datetime.now().strftime("%y%m%d-%H%M%S") + ".txt"
header = "PDB_file, probe, C666SG_dist, C666_access, C666_clashes, all_clashes, fract_COM, BD_L2, AC_L2, del_dimer_COM\n".replace(", ","\t")
l = open(logfile, 'w+')
l.write(header)
l.close()

#fixed parameters
# enter the list of pdb files here
source_PDBs = ["LBD_tet", "fw_noh", "loose", "tight"]
# use the right residue number for each input PDB for Cys and Linker site (666 is 154 in LBD-only construct)
residues = [666, 154, 154, 154]
gates = [632, 120, 120, 120]
dimer_one_subs = ["A", "D"]
dimer_two_subs = ["C", "B"]

probe = 3           #the probe size for areaimol (bigger than default - M1M in mind
minClose = 2.2      #the minimum distance for a clash


for seed, residue, gate in zip(source_PDBs, residues, gates):
    pymol.cmd.reinitialize()
    seedPBD = seed + ".pdb"
    pymol.cmd.load(seedPBD)

    print "Measurements for {0}.pdb".format(seed)
    results = []

    #Seed is the stem of the filename and therefore also the automatic object name in PyMOL.
    C666SG_dist = pymol.cmd.distance("/" + seed + "//A/{0}/SG".format(residue), "/" + seed + "//C/{0}/SG".format(residue))
    
    BD_L2 = pymol.cmd.distance("/{0}//B/{1}/CA".format(seed, gate), "/{0}//D/{1}/CA".format(seed, gate))
    AC_L2 = pymol.cmd.distance("/{0}//A/{1}/CA".format(seed, gate), "/{0}//C/{1}/CA".format(seed, gate))

    d1_com = pymol.cmd.centerofmass("Chain " + "+".join(dimer_one_subs))
    d2_com = pymol.cmd.centerofmass("Chain " + "+".join(dimer_two_subs))
    dist = [(a - b) ** 2 for a, b in zip(d1_com, d2_com)]
    del_dimer_COM = math.sqrt(sum(dist))
    fract_COM = 0   #not relevant, no movement

    all_clashes = getClashes (seedPBD, "NCout", "A,D", "B,C", minClose)
    C666_clashes = getClashes (seedPBD, "NC_resi_out", "A/{0}/SG".format(residue), "B,C", probe)
    C666_access = getAccess (seedPBD, "AIout", 4, "CYS A {0}".format(residue))

        
    measurements = [probe, C666SG_dist, C666_access, C666_clashes, all_clashes, fract_COM, BD_L2, AC_L2, del_dimer_COM]
    hS = hitCalculation(measurements) #calculate but don't store
    model_details = "C666 SG dist. (Ang.), accessible area (to probe of {m[0]} Ang): {m[1]:.1f}, {m[2]:.2f}\n\
    Clashes (C666, all): {m[3]}, {m[4]}\n\
    Dimer Centre of mass distance change, fractional change : {m[8]:.1f}, {m[5]:.1f}\n\
    Distances (Ang) for BD L2, AC L2:  {m[6]:.1f}, {m[7]:.1f}\n\
    Hitscore = {h}\n".format(m=measurements, h=hS)
            
    print model_details
    results.append(measurements)

    l = open(logfile, 'a+')
    l.write(seed + ".pdb\t" + formatLine(measurements))
    l.close()
     
