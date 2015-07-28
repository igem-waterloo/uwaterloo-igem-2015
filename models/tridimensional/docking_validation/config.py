import os


DIR_BASE = "~/working/"  # this should be changed on windows systems

# ====================================================================== #
# Directory layout                                                       #
# ====================================================================== #
#      /--- PDB ---/   pdb files (protein-dna structures) go here        #
# BASE /- RESULTS -/   results and simulation scores go here             #
#      /- SCRIPTS -/   all scripts go here (should be flat)              #
# ====================================================================== #

DIR_PDB = os.path.join(DIR_BASE, "pdb")
DIR_RESULTS = os.path.join(DIR_BASE, "results")
DIR_SCRIPTS = os.path.join(DIR_BASE, "scripts")


# ====================================================================== #
# PDB:                                                                   #
# ====================================================================== #
#                                                                        #
#     /- special -/    any special pdb files, esp. the progenitor pdb    #
#     |                                                                  #
# pdb /- mutants -/-- mutant_0 --/-- template --/   1 pdb                #
#                 |              /-- variants --/   256 pdbs             #
#                 |                                                      #    
#                 /-- mutant_1 --/-- template --/   1 pdb                #
#                 |              /-- variants --/   256 pdbs             #
#                 |                                                      #
#                 etc.                                                   #
#                                                                        #
# Every mutant gets its owned subdirectory of pdb files based on index.  #
# For example, if index=12, then the mutant has a template associated    #
# with applying the index=12 mutations from our mutant database to       #
# the original pdb file in  "pdb/special/".                              #
# This template pdb is left in:  "pdb/mutants/mutant_12/template/"       #
# Once the template is created, 64 or 256 pam variants are created and   #
# stored in:  "pdb/mutants/mutant_12/variants/".                         #
# ====================================================================== #

DIR_PDB_MUTANTS = os.path.join(DIR_BASE, "mutants")
DIR_PDB_SPECIAL = os.path.join(DIR_BASE, "special")
PATH_PDB_SPECIAL_ORIGINAL = os.path.join(DIR_PDB_SPECIAL, "4UN3.original.pdb")
PATH_PDB_SPECIAL_VQR = os.path.join(DIR_PDB_SPECIAL, "4UN3.mutant_VQR.pdb")
PATH_PDB_SPECIAL_EQR = os.path.join(DIR_PDB_SPECIAL, "4UN3.mutant_EQR.pdb")
FOLDERNAME_MUTANT_TEMPLATE = "template"
FOLDERNAME_MUTANT_VARIANTS = "variants"
MUTANT_TEMPLATE_PREFIX = "4UN3.mutant_"


# ====================================================================== #
# RESULTS:                                                               #
# ====================================================================== #
#                                                                        #
#         /- batch -/   validation runs for specific tests/analysis      #
#         |                                                              #
# results /- mutant_scores -/- scores_0 -/- run_0 -/   256 txt, 1 csv    #
#                           |            /- run_1 -/   256 txt, 1 csv    #
#                           |                ...                         #
#                           |            /- run_k -/   256 txt, 1 csv    #
#                           |                                            #
#                           /- scores_1 -/- run_0 -/   256 txt, 1 csv    #
#                           |            /- run_1 -/   256 txt, 1 csv    #
#                           |                ...                         #
#                           |            /- run_k -/   256 txt, 1 csv    #
#                           |                                            #
#                           etc.                                         #
#                                                                        #
# The results directory mirrors the pdb directory in that each mutant    #
# gets its own folder: "results/mutant_scores/scores_i" where i is the   #
# mutant index. Within scores_i are k "run_j" folders. We repeat runs    #
# to characterize the stochastic elements of the monte carlo simulation. #
# ====================================================================== #

DIR_RESULTS_VALIDATION = os.path.join(DIR_RESULTS, "batch")
DIR_RESULTS_MUTANT_SCORES = os.path.join(DIR_RESULTS, "mutant_scores")


# ====================================================================== #
# SCRIPTS:                                                               #
# ====================================================================== #
#                                                                        #
# scripts /   all scripts here                                           #
#                                                                        #
# This is currently a flat directory that may have some structure later. #
# ====================================================================== #

MUTANT_DATABASE_PATH = os.path.join(DIR_SCRIPTS, "mutant_database.py")
