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
##
# to get access to CCP4 the script should be run from the Terminal with CCP4 sourced.
# use NCONT to count close contacts for this residue ("buried")
# use AREAIMOL to see if residue has been buried
# [proposed] use CNS to do simulated annealing on models and tidy them up (like in Morph)
# <minimize powell> is the command (from morph_dist.inp)
##
# Tell PyMOL to launch quiet (-q), and with no GUI (-c)
# Make sure PyMOL has finished to launch before using any PyMOL modules.

# Beware - many of the numbers here may have to be tweaked for your application
# for example, Cys max accessible area is 30 ang^2
# would be different for other residues
# different structures may require other convergence parameters

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
source = "LBD_tet"
residue_bury = 666
perp_axis = [0,     1,     0] #y??
para_axis = [1,     0,     0] #x??
dimer_one_subs = ["A", "D"]
dimer_two_subs = ["C", "B"]
n = 5               #number of new structures to generate for each search
probe = 3           #the probe size for areaimol (bigger than default - M1M in mind
minClose = 2.2      #the minimum distance for a clash
minBD = 55          # need to maintain BD distance for channel activation
minAC = 20          # need to maintain AC distance to keep dimers in plane

#loop and search for 50 models -takes about 4 hours (180305)
for mod in range(50):
    #reset adjustable params for a fresh search
    best_file = None
    best_hit = 30       #max expousre of Cys 666
    search = True
    round = 0           #begin at zeroth round of search
    seed = source       #reset to use original file
    rot_range = 20      #degrees
    tra_range = 10      #Angstroms
    lost = 0            #not lost in current search yet
    no_improve = False  #improvement is possible
    
    while search:
        pymol.cmd.load(seed + ".pdb")
        round += 1
        print "Round {0}: best hit score {1}".format(round , best_hit)
        results = []
        for trial in range(n):
            trial_file = "trial" + str(trial)
            trial_PDBname = trial_file + ".pdb"
            print "Trial {0}".format(trial)
            #pymol is reinitialized each call
            seedCOM, rotPL, rotPP, tx, ty = rearrange_subunits (seed, trial_file, dimer_one_subs, dimer_two_subs, rot_range, tra_range, para_axis, perp_axis)

            #Seed is the stem of the filename and therefore also the automatic object name in PyMOL.
            C666SG_dist = pymol.cmd.distance("(/" + seed + "//A/666/SG)", "(/" + seed + "//C/666/SG)")
            BD_L2 = pymol.cmd.distance("(/" + seed + "//B/632/CA)", "(/" + seed + "//D/632/CA)")
            AC_L2 = pymol.cmd.distance("(/" + seed + "//A/632/CA)", "(/" + seed + "//C/632/CA)")

            d1_com = pymol.cmd.centerofmass("Chain " + "+".join(dimer_one_subs))
            d2_com = pymol.cmd.centerofmass("Chain " + "+".join(dimer_two_subs))
            dist = [(a - b) ** 2 for a, b in zip(d1_com, d2_com)]
            del_dimer_COM = math.sqrt(sum(dist))
            fract_COM = del_dimer_COM / seedCOM

            all_clashes = getClashes (trial_PDBname, "NCout", "A,D", "B,C", minClose)
            C666_clashes = getClashes (trial_PDBname, "NC_resi_out", "A/666/SG", "B,C", probe)
            C666_access = 30.0 #default high value. When it was only 1, some bad structures got accepted
            
            #AREAIMOL is slow
            #no obvious way to reduce the input/output to speed up
            #so only launch it if there are not clashes
            if all_clashes < 10 and C666_clashes == 0:
                print "Few enough clashes (all: {0}, resi {2}: {1}) so use AREAIMOL to see if {2} buried".format(all_clashes, C666_clashes, residue_bury)
                C666_access = getAccess (trial_PDBname, "AIout", 4, "CYS A 666")
                if C666_access > 30:
                    C666_access = 30    #full access
        
            measurements = [probe, C666SG_dist, C666_access, C666_clashes, all_clashes, fract_COM, BD_L2, AC_L2, del_dimer_COM]
            hS = hitCalculation(measurements) #calculate but don't store
            model_details = "Movements- rotPL {0:.1f} rotPP {1:.1f} tx {2:.1f} ty {3:.1f}\n\
    C666 SG dist. (Ang.), accessible area (to probe of {m[0]} Ang): {m[1]:.1f}, {m[2]:.2f}\n\
    Clashes (C666, all): {m[3]}, {m[4]}\n\
    Dimer Centre of mass distance change, fractional change : {m[8]:.1f}, {m[5]:.1f}\n\
    Distances (Ang) for BD L2, AC L2:  {m[6]:.1f}, {m[7]:.1f}\n\
    Hitscore = {h}\n".format(rotPL, rotPP, tx, ty, m=measurements, h=hS)
            
            print model_details
            results.append(measurements)

        ### work out if we have found a good model
        hitScoreAvg = 0
        for index, measure in enumerate(results):
            
            hit = hitCalculation(measure)   #simple so recalculate
            hitScoreAvg += float(hit) / n
            #show model statistics
            print "{3}: hitScore {0} (best {1}), measures - {2}".format(hit, best_hit, formatLine(measure), index)
            if measure[6] > minBD and measure[7] > minAC:
                if (hit < best_hit and measure[5] < 1.02) or (hit == best_hit and measure [5] < 1):
                    #dimers should approach or at worse move slightly apart and stay roughly in plane...
                    seed = "trial" + str(index)     #use new best model as the seed for subsequent search
                    print "{} is new best model!".format(seed)
                    best_hit = hit
                    best_measure = measure  #store the measurements to write out later
                    
                    # reduce the step size to more effectively search around best point
                    rot_range /= 1.4
                    tra_range /= 1.4
                    print "Step sizes reduced to {0:.2f} degrees rotation and {1:.2f} angstrom translation".format(rot_range, tra_range)
                    best_file = "best" + str(round)
                    os.system("cp {0} {1}.pdb".format(seed + ".pdb", best_file)) # keep the best fit so far
            else:
                print "Bad BD or AC distance"

        #envisage that most of the movements generate a lot of clashes - dimers must be close.
        if hitScoreAvg > 300:
            rot_range /= 1.4
            tra_range /= 1.4
            print "Too many clashes, step sizes reduced to {0} degrees rotation and {1} angstrom translations".format(rot_range, tra_range)

        #if no clashes occur for several rounds, it's likely that the dimers are too far apart
        #therefore: follow COM reduction; don't improve best model, but alter seed
        if hitScoreAvg == 30:
            #the average score for 5 completely free models (max CYS access at 666 set at 30 Ang^2)
            lost += 1
            print "Are dimers lost in space?, {}/3 rounds before pushing back together"
        else:
            lost = 0    #reset because there was a clash this time
        
        if lost > 3:
            print "Lost? Trying to move in direction of dimer center-of-mass displacement reduction"
            #move in direction of COM reduction
            best_COM_red_index = None
            best_COM_reduction = 1
            for index, measure in enumerate(results):
                if measure[5] < best_COM_reduction:
                    best_COM_reduction = measure[5]
                    best_COM_red_index = index

            if best_COM_red_index:
                seed = "trial" + str(best_COM_red_index)
                print "changed seed to move dimers back together {} {}".format(seed, best_COM_reduction)

        #catch each reason the model cannot be betttered separately and produce message
        if tra_range < 0.3:
            no_improvement_message = "converged, translation < 0.3A"
            no_improve = True
        if round == 30:
            no_improvement_message = "30 rounds passed without converging - give up"
            no_improve = True
            if best_file == None:
                no_improvement_message += ", sadly no best model found."
        if best_hit == 0.0:
            no_improvement_message = "perfectly buried: best_hit = 0"
            no_improve = True
         
        #obvious that no improvement can occur - write out best model or give up
        if no_improve == True:
            print no_improvement_message
            search = False
            #it could be that we reach round 30 without finding a good file
            if best_file != None:
                timestamp = datetime.now().strftime("%y%m%d-%H%M%S")
                #write out best with time-stamp and repeat
                os.system("mv {0}.pdb {1}.pdb".format(best_file, "best"+timestamp))
                # there must be a best measure - this will write out horribly - format and add data
                # tidy writing and closing to see progress...
                l = open(logfile, 'a+')
                l.write("best" + timestamp + ".pdb\t" + formatLine(best_measure))
                l.close()

    print "Search {} ended, hitScore: {}".format(mod, best_hit)     #in case there is no best file, don't try


