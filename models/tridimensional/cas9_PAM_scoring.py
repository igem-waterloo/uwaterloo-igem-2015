"""
IMPORTANT
Modify script to point to the location of iGEM Dropbox on your computer. Thanks! :)

- Kevin

"""


from time import time
import os
from rosetta import *

# Navigate to where the PDBs are located
os.chdir("/home/aruu/Documents/Dropbox/Waterloo iGEM 2015") # Change to location of your iGEM Dropbox
os.chdir("Math Modelling/cas9_modification")

# Initialize PyRosetta
init()

# Variations to account for
dna_nts = "acgt"
programs = ["Chimera", "3DNA"]


# If we wanna parallelize+distribute this, probably easiest just for each person to
# hard code what range they are responsible for (to avoid everyone processing the
# same 4 tasks, no way to tell what has already been started by someone else)

for i in range(0, 64):
	# Creates a three nucleotide PAM sequence, aaa then, aac, aag, etc.
	variant = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]

	for program in programs:
		results_file = "results/" + variant + "_" + program + ".txt"
		# Next if results already exist
		if os.path.isfile(results_file): continue

		print "Running for: " + variant + "_" + program
		# Start keeping track of runtime
		time_init = time()
		test_pose = pose_from_pdb(program + "/4UN3." + variant + ".pdb")

		# Create a custom FoldTree and replace the one which our Pose is using
		# TODO: Find out what exactly FoldTrees do, are there better ones for our purpose
		pose_fold_tree = FoldTree(test_pose.total_residue())
		pose_fold_tree.new_jump(1, 1307, 1306)
		test_pose.fold_tree(pose_fold_tree)

		# Full-atom scoring function
		fa_score = get_fa_scorefxn()

		# DNA scoring function (from demo/D110_DNA_interface.py)
		dna_score = create_score_function('dna')
		dna_score.set_weight(fa_elec , 1)

		# Docking protocol
		# TODO: Look into how to modify this Docking protocol
		docking = DockMCMProtocol()
		docking.set_scorefxn(dna_score)
		docking.set_partners("B_ACD")


		fa_score_init = fa_score(test_pose)
		dna_score_init = dna_score(test_pose)
		docking.apply(test_pose)
		fa_score_final = fa_score(test_pose)
		dna_score_final = dna_score(test_pose)
		time_final = time()

		print "DONE RUNNING"
		f = open(results_file, 'a')		# Append results in case file exists
		text = "Final FA Score:  %8.3f\n" % (fa_score_init)
		f.write(text)
		text = "Init DNA Score:  %8.3f\n" % (dna_score_init)
		f.write(text)
		text = "Final FA score:  %8.3f\n" % (fa_score_final)
		f.write(text)
		text = "Final DNA score: %8.3f\n" % (dna_score_final)
		f.write(text)
		text = "Time taken:      %8.3f\n" % (time_final- time_init)
		f.write(text)
		f.write("\n")
		f.close()
		print "DONE WRITING"