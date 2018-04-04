#!/usr/bin/env python

import os
import pymol
import random
from pymol import cmd
import math

def formatLine (list):
    #get a tab delimited string
    out = ""
    for e in list:
        if type(e) is not float:
            out += str(e) + '\t'
        else:
            out += "{:.1f}".format(e) + '\t'
    return out + "\n"

def getClashes (inputPDB, outfile, src, trg, min):
    #find clashes using NCONT from CCP4 - must be sourced in shell (.profile on Mac)
    NC_outfile = outfile+".txt"
    
    cmd  = """ncont xyzin {0} > {1} <<eof
        source  {2}
        target  {3}
        maxdist {4}
        end
        """.format(inputPDB, NC_outfile, src, trg, min)

    # run NCONT
    os.system(cmd)
    
    # parse NCONT output file
    with open(NC_outfile, "r") as g:
        searchlines = g.readlines()
        
        #start by assuming no clashes
        clashes = 0
        for j, line in enumerate(searchlines):
            
            if "contacts found:" in line:
                #print j, line
                #the third element is the one with the number of contacts
                #with older version of NCONT, element 2 was the right one
                clashes = int(line.split(' ')[2])

    return clashes

def getAccess (inputPDB, outfile, probe_radius, residue):
    #find residues that can't be accessed by a large probe using AREAIMOL from CCP4 - must be sourced in shell (.profile on Mac)
    #residue : a string that specifies resi and chain like "CYS A 666"
    #
    
    AI_outfile = outfile + ".txt"
    cmd = """areaimol XYZIN {0} > {1} <<eof
        VERB      ! Verbose output
        PROBE {2}   ! chose probe radius
        END
        eof""".format(inputPDB, AI_outfile, probe_radius)
    os.system(cmd)

    with open(AI_outfile, "r") as g:
        #start by assuming full access
        access = 30
        searching = True
        while searching:
            for line in g:
            
                left, sep, right = line.partition(residue)
                if sep:             #only true for line where [residue] is found
                    #print line
                    #accessArea =float(right[8:16]) #normalized GXG value
                    accessArea =float(right[:8]) #full area
                    print "Access area", accessArea, line
                    searching = False   #we want the first instance only.
                    break
                    
        return accessArea

def hitCalculation(m):
    #m :
    #probe, C666SG_dist, C666_access, C666_clashes, all_clashes, fract del dimer com, BD_L2, AC_L2, del_dimer_COM
    #low score is better - perfect score of zero would be no access to a non-clashing C666
    #score is sum of C666 exposed absolute area, C666 clashes and the total # of clashes
    score = float(m[2])  + float(m[3]) + float(m[4])
    return score

def rearrange_subunits(seed, outpdb, group1, group2, rot_range, tra_range, para_axis, perp_axis):
    # group1 and group2 should be a list of chains that should be moved.
    
    # reinitialize
    pymol.cmd.reinitialize()
    #take seed
    pymol.cmd.load(seed + ".pdb")
    
    sgroup1 = "+".join(group1)
    sgroup2 = "+".join(group2)
    
    pymol.cmd.orient(seed)      #with very badly dispersed models, this fails
    
    pymol.cmd.select("g1", "chain " + sgroup1)
    pymol.cmd.select("g2", "chain " + sgroup2)
    
    g1_com = pymol.cmd.centerofmass("g1")
    g2_com = pymol.cmd.centerofmass("g2")
    dist = [(a - b) ** 2 for a, b in zip(g1_com, g2_com)]
    seedCOM = math.sqrt(sum(dist))
    
    pp = random.uniform(-rot_range, rot_range)
    pl = random.uniform(-rot_range, rot_range)
    
    tx = random.uniform(-tra_range, tra_range)
    ty = random.uniform(-tra_range, tra_range)
    
    #rotate each dimer in turn
    pymol.cmd.rotate(perp_axis, pp, "g1")
    pymol.cmd.rotate(perp_axis, -pp, "g2")
    
    pymol.cmd.rotate(para_axis, pl, "g1")
    pymol.cmd.rotate(para_axis, -pl, "g2")
    
    # translate symmetrically in xy plane a bit
    pymol.cmd.translate("[{0}, {1} , 0]".format(tx, ty), "g1")
    pymol.cmd.translate("[{0}, {1} , 0]".format(-tx, -ty), "g2")
    
    pymol.cmd.save(outpdb + ".pdb", seed)         #these will get overwritten each round

    return seedCOM, pp, pl, tx, ty


