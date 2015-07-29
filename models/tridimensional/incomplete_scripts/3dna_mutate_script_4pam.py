import os
# import subprocess to run cli
import subprocess

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
for i in listOfBases:
    for j in listOfBases:
        for k in listOfBases:
            for l in listOfBases:
                PAM = i+j+k+l
                # TGG PAM site already exists, let's skip it
                mutateNucleotide(fileNameOriginal,PAM)
                print "mutate PAM to " + PAM
                # Save and close all files
