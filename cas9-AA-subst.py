#///////////////////////**READ ME FIRST**//////////////////////////////////
# Point of this script is to substitute amino acids into a given pose 
#  protein structure in pyrosetta.
# Please change dropbox folder location accordingly + make sure that the 
#  script's current location is the root of PyRosetta (need to import stuff)
#//////////////////////////////////////////////////////////////////////////
## Raw file date modified: July 15th/2014

## importing modules required for script
import os
from rosetta import *
## difference between resfiles and mutate_residue is that the former takes into account 
##  rotamers (aka conformers) whereas the latter doesn't (for our case mutate_residue is sufficient)
##  refer to link for further documentation and/or info:
##  http://graylab.jhu.edu/pyrosetta/downloads/documentation/Workshop6_PyRosetta_Packing_Design.pdf
from toolbox import generate_resfile_from_pdb # generate mutations using resfiles
from toolbox import mutate_residue # generate mutations using mutate_residue  

## changing directory to where PDB's are located (aka where PDB files are located )
os.chdir("~\Dropbox\Waterloo-iGEM-2015") #alter to your specific dropbox path
os.chdir("\Math Modelling\cas9_modification") ##where the WT cas9 should be located???
## not sure if completely correct???? add changes if not.

## initializing rosetta:
rosetta.init()

# import cleaning module for PDB to be usable
from toolbox import cleanATOM

cleanATOM("\4UN3.pdb") # cleaned PDB file to use for analysis
var_pose = pose_from_pdb("\4UN3.pdb") # initial pose created from clean pdb

#inputted residue number of interest
Num = raw_input("enter residue number:\n")

for i in range(0, 20):
	# list of Amino Acids to substitute
	AA_lst = "ACDEFGHIKLMNPQRSTVWY"
	AA_var = AA_lst[i]

	var_pose = pose_from_pdb("4UN3." + AA_var + ".clean.pdb")
	mutate_residue(var_pose, Num , AA_var) # where Num = residue number to substitute AA

	# for sanity checking purposes: prints out changed 4UN3 pdb pose protein profile 
	#  (sequence, # of res, what is located at Num residue - if the substitution occured)
	print var_pose
	print var_pose.sequence()
	print "Protein has", var_pose.total_residue(), "residues."
	print var_pose.residue(Num).name() # where Num is a residue number


