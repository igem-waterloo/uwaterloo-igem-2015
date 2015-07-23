"""
IMPORTANT
Modify script to point to the location of iGEM Dropbox on your computer. Thanks! :)
- Kevin
"""

from time import time
import os
from multiprocessing import Pool
# Navigate to where the PDBs are located
os.chdir("C:\Users\mark\My Documents\Dropbox\Mark - iGEM\Cas9_rational_design")

# Set up pyrosetta
from rosetta import *
# option reduces amount of output, nice for long scripts
rosetta.init(extra_options="-mute all")
# Variations to account for
dna_nts = "acgt"
programs = ["Chimera", "3DNA"]
""" It's only running a single example, but I think I found an easier way to 
paralellize/divide up the scripts. This is list of numbers now, but it would
be simple to adjust to using the actual PAM sequence
- Mark"""
listOfPAMsToRun = range(0,64)
threads = 2

def scorer(i):
    # Creates a three nucleotide PAM sequence, aaa then, aac, aag, etc.
    variant = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]
    for program in programs:
        # Start keeping track of runtime
        time_init = time()
        results_file = 'results_new.txt'
        print "Running for: " + variant + "_" + program
        original_pose = pose_from_pdb(program + "\\4UN3." + variant + ".pdb")
        setup_foldtree(original_pose, 'B_CD', Vector1([1]))

        test_pose = Pose()
        test_pose.assign(original_pose)

        # DNA scoring function (from demo/D110_DNA_interface.py)
        dna_score = create_score_function('dna')
        dna_score.set_weight(fa_elec , 1)

        # Full-atom scoring function
        fa_score = get_fa_scorefxn()

        # Docking protocol
        docking = DockMCMProtocol()
        docking.set_scorefxn(dna_score)
        docking.set_scorefxn_pack(fa_score)
        docking.set_partners("B_ACD")

        dna_score_init = dna_score(test_pose)
        fa_score_init  = fa_score(test_pose)
        docking.apply(test_pose)
        dna_score_docking = dna_score(test_pose)
        fa_score_docking  = fa_score(test_pose)
        time_final = time()
        print "Finished Docking Variant: " + variant + "_" + program
        # 'a+' appends results if file exists. Output formatted for single csv file
        f = open(results_file, 'a+')    
        text = variant[0] + "," + variant[1] + "," + variant[2] + "," + program + ","
        f.write(text)
        text = "%8.3f," % (dna_score_init)
        f.write(text)
        text = "%8.3f," % (dna_score_docking)
        f.write(text)
        text = "%8.3f," % (fa_score_init)
        f.write(text)
        text = "%8.3f," % (fa_score_docking)
        f.write(text)
        text = "%8.3f,\n" % (time_final- time_init)
        f.write(text)
        f.close()
        print "Finished Writing Score for Variant: " + variant + "_" + program

if __name__ == '__main__':
    # Uncomment the next two lines to parallize, can just pass any list to 
    p = Pool(processes=threads)
    p.map(scorer,listOfPAMsToRun)