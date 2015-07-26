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

# Variations to account for
dna_nts = "acgt"
programs = ["Chimera", "3DNA"]

threads = 2
""" It's only running a single example, but I think I found an easier way to 
paralellize/divide up the scripts. This is list of numbers now, but it would
be simple to adjust to using the actual PAM sequence
- Mark"""

listOfPAMsToRun = [0,5,10,19,20,58]
jobs = 1

def scorer(i):
    # Creates a three nucleotide PAM sequence, aaa then, aac, aag, etc.
    variant = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]
    rosetta.init(extra_options="-mute all")
    for program in programs:
        results_file = 'results.txt'
        print "Running for: " + variant + "_" + program
        # Start keeping track of runtime
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

        job_output = 'results\\' + program + '_' + variant

        jd = PyJobDistributor( job_output, jobs, dna_score)
        counter = 0
        while not jd.job_complete:
            test_pose.assign(original_pose)
            counter +=1
            test_pose.pdb_info().name(job_output + '_' + str(counter))
            time_init = time()
            dna_score_init = dna_score(test_pose)
            docking.apply(test_pose)
            dna_score_docking = dna_score(test_pose)
            time_final = time()
            test_pose.pdb_info().name(job_output + '_' + str(counter) + '_fa')
            # !!!!!IMPORTANT!!!!!
            # Do not comment out the line below, or the script gets stuck in a 
            # loop scoring the same pdb repeatedly, regardless of the job number
            jd.output_decoy(test_pose)
            print "DONE RUNNING"
            f = open(results_file, 'a+')     # Append results in case file exists
            text = variant[0] + "," + variant[1] + "," + variant[2] + "," + program + ","
            f.write(text)
            text = "%8.3f," % (dna_score_init)
            f.write(text)
            text = "%8.3f," % (dna_score_docking)
            f.write(text)
            text = "%8.3f,\n" % (time_final- time_init)
            f.write(text)
            f.close()
            print "DONE WRITING"

if __name__ == '__main__':
    # Uncomment the next two lines to parallize, can just pass any list to 
    p = Pool(processes=threads)
    p.map(scorer,listOfPAMsToRun)
