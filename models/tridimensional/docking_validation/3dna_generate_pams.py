from __future__ import print_function

import sys
import os
import argparse
import math
import errno

import subprocess
from constants import PAM_TEMPLATE_SEQUENCE, DNA_ALPHABET

# Nucleotides we want to mutate are located at positions 5,6,7 in chain D (PAM,NGG)
# and 8,7,6 in chain C (target strand, NCC).
positionIndex = { 0 : "5", 1 : "6", 2 : "7", 3 : "8"}
# Create a dictionary mapping corresponding positions to each other.
positionPairs = { "5" : "8", "6" : "7", "7" : "6", "8" : "5"}

# Create dictionary mapping valid base pairs to each other.
basePairs = {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'}

# Function to mutate base pairs
def mutateNucleotide(fileNameOriginal,PAM):
    mutationString = createMutationString(PAM)
    print mutationString
    subprocess.call(["mutate_bases", mutationString, fileNameOriginal, "4UN3." + PAM + ".pdb"])

# Function to determine command line output for mutations and positions for
# call call of mutate_base
def createMutationString(PAM):
    # Set up some variables
    mutationString = ""
    # Loop through the PAM sequence and find positions to mutate
    for n in range(len(PAM)):
        pos = positionIndex[n]
        base = PAM[n]
        complementaryBase = basePairs[base]
        complementaryPosition = positionPairs[pos]
        # If the first base is T, don't mutate it
        if base== 't' and n == 0:
            continue
        # If the 2nd or 3rd bases are G, don't mutate them
        if all([base == 'g', any([n == 1, n == 2])]):
            continue
        mutationString = mutationString + "c=d s=" + pos + " m=D" + base + \
                "; c=c s=" + complementaryPosition + " m=D" + complementaryBase + "; "
    mutationString = mutationString[0:-2]
    return mutationString
# Going to assume file"4UN3.clean.pdb" is always used and in the same directory.
# This file has a TGG PAM sequence already present in it
fileNameOriginal = "4UN3.original.pdb"

# Iterate through all possible PAM sites
listOfBases = ['a','c','g','t']

if __name__ == '__main__':
    # create parser and parse arguments
    parser = argparse.ArgumentParser(description='Generate PDBs with new PAM sites based on input PDB with Cas9 variant')
    parser.add_argument('-n', '--num_pams', metavar='N', type=str, # string type because of how this script must be run by Chimera
                        help='how many PAMs in total to run. 64 = all PAMs of length 3, 256 = all PAMs of length 4')
    parser.add_argument('-i', '--input_pdb', metavar='F', type=str,
                        help='input PDB file containing Cas9 mutant of interest')
    parser.add_argument('-o', '--output_dir', metavar='D', type=str,
                        help='path to output directory for new PDBs')
    args = parser.parse_args()

    assert args.num_pams is not None
    assert args.input_pdb is not None
    assert args.output_dir is not None
    args.num_pams = int(args.num_pams) # convert string from command line to int
    pam_length = int(math.log(args.num_pams, 4))
    assert 64 == args.num_pams or 256 == args.num_pams
    assert os.path.isfile(args.input_pdb)

    try: # check existence again to handle concurrency problems
        os.makedirs(args.output_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(args.output_dir):
            pass
        else: raise


    for i in listOfBases:
        for j in listOfBases:
            for k in listOfBases:
                for l in listOfBases:
                    PAM = i+j+k+l
                    # TGG PAM site already exists, let's skip it
                    mutateNucleotide(fileNameOriginal,PAM)
                    print "mutate PAM to " + PAM
                    # Save and close all files
