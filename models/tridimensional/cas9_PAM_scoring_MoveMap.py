"""
IMPORTANT
Modify script to point to the location of iGEM Dropbox on your computer. Thanks! :)
- Kevin
"""

from time import time
import os
from rosetta import *
from multiprocessing import Pool
# Navigate to where the PDBs are located
os.chdir("C:\Users\mark\My Documents\Dropbox\Mark - iGEM") # Change to location of your iGEM Dropbox
os.chdir("Cas9_rational_design")

# Initialize PyRosetta
init()

# Variations to account for
dna_nts = "acgt"
programs = ["Chimera", "3DNA"]

threads = 2
""" It's only running a single example, but I think I found an easier way to 
paralellize/divide up the scripts. This is list of numbers now, but it would
be simple to adjust to using the actual PAM sequence
- Mark"""

listOfPAMsToRun = [0,5,10,19,20,58]

def scorer(i):
    # Creates a three nucleotide PAM sequence, aaa then, aac, aag, etc.
    variant = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]

    for program in programs:
        results_file = "results.txt"

        print "Running for: " + variant + "_" + program
        # Start keeping track of runtime
        time_init = time()
        original_pose = pose_from_pdb(program + "\\4UN3." + variant + ".pdb")
        test_pose = Pose()
        test_pose.assign(original_pose)

        # Create a movemap specify what AAs to move
        pose_move_map = MoveMap()
        # Set all residues to immutable to start
        pose_move_map.set_bb(False)
        pose_move_map.set_chi(False)
        # This should say 1117-1387 for the PI domain, changed for testing
        # Wanted to make sure it was within the Protein strand
        pose_move_map.set_bb_true_range(1117,1287)
        pose_move_map.set_chi_true_range(1117,1287)

        # DNA scoring function (from demo/D110_DNA_interface.py)
        dna_score = create_score_function('dna')
        dna_score.set_weight(fa_elec , 1)

        # Full-atom scoring function
        fa_score = get_fa_scorefxn()
        """
        Ignore this, I'm trying to figure out how a bunch of the different
        Packers work
        # Relax the atom
        fa_score_init = fa_score(test_pose)
        relax = ClassicRelax()
        relax.set_scorefxn(fa_score)
        relax.set_movemap(pose_move_map)
        relax.apply(test_pose)
        fa_score_relax = fa_score(test_pose)
        """ 
        # Docking protocol
        # TODO: Look into how to modify this Docking protocol
        docking = DockMCMProtocol()
        print "Loaded Docking"
        docking.set_scorefxn(dna_score)
        print "set docking score function"
        docking.set_partners("B_ACD")
        print "Set Partners"
        docking.set_move_map(pose_move_map)
        # This -should- restrict repacking of residues to specific range, but
        # doesn't seem to work
        print "Set Move Map"
        # Initial Scoring
        dna_score_init = dna_score(test_pose)
        print "Done Initial Scoring"
        # Score after docking
        print "Begin Docking"
        # The script is currently crashing in here, I can't find out why
        docking.apply(test_pose)
        print "Finished Docking"
        dna_score_docking = dna_score(test_pose)

        time_final = time()

        print "DONE RUNNING"
        f = open(results_file, 'a+')     # Append results in case file exists
        text = variant[0] + "," + variant[1] + "," + variant[2] + "," + program + ","
        f.write(text)
        text = "%8.3f," % (dna_score_init)
        f.write(text)
        text = "%8.3f," % (dna_score_docking)
        f.write(text)
        text = "%8.3f," % (time_final- time_init)
        f.write(text)
        f.close()
        print "DONE WRITING"

if __name__ == '__main__':
    # Uncomment the next two lines to parallize, can just pass any list to 
    #p = Pool(processes=threads)
    #p.map(scorer,listOfPAMsToRun)
    scorer(0)