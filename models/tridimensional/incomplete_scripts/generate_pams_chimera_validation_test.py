import argparse
import errno
import math
import os
import sys

import chimera
from chimera import runCommand
from chimera import replyobj
HOME = os.getenv("HOME")
sys.path.append(HOME+'/github/uwaterloo-igem-2015/models/tridimensional/docking_validation')
from constants import PAM_TEMPLATE_SEQUENCE, DNA_ALPHABET
from utility import pam_string_from_int 
'''
Mark Lubberts - 01/08/2015
I'm checking if the higher scoring in the tgg variants is due to the fact that
they aren't mutated, rather than actually binding better. tgg bases are switched
with identical nucleotide instead of being skipped, to ensure that scoring is 
nucleotide specific.
'''

def mutate_nt(pam_idx, base):
    """Mutate a base pair at pam_idx to base
    Args:
        pam_idx: the index of the PAM nucleotide to be modified
        base: what the nt at pam_idx will be changed to
    """
    position_idx = { 0 : "5.d", 1 : "6.d", 2 : "7.d", 3 : "8.d"} # Nucleotides we want to mutate are located at positions 5,6,7 in chain D (PAM,NGG) and 8,7,6 in chain C (target strand, NCC).
    position_pairs = { "5.d" : "8.c", "6.d" : "7.c", "7.d" : "6.c", "8.d" : "5.c"} # Create a dictionary mapping corresponding positions to each other.
    base_pairs = {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'} # Create dictionary mapping valid base pairs to each other.

    pos = position_idx[pam_idx]
    complement_base = base_pairs[base]
    complement_pos = position_pairs[pos]
    runCommand("swapna " + base + " : " + pos )
    runCommand("swapna " + complement_base + " : " + complement_pos)


def generate_pam_variant_chimera(pam, input_pdb, output_dir):
    """Generate a new PDB file
    Args:
        pam: new PAM sequence to be listed in name
        template_pdb_basename: original fine basename
        output_dir: output directory
    """
    template_pdb_basename = os.path.basename(input_pdb)
    new_pdb_path = os.path.join(output_dir, template_pdb_basename[:-4] + "." + pam + ".pdb")
    runCommand("write 0 " + new_pdb_path)
    runCommand("close all")


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

    for i in xrange(args.num_pams):
        pam = pam_string_from_int(i, pam_length)

        # open up the file again each time, for now
        runCommand("open " + args.input_pdb)

        # loop through the PAM sequence and mutate positions
        for pam_idx in xrange(pam_length):
            '''This has been commented out so that all nucleotides are changed
            regardless of whether the match the original pdb or not.'''
            # If the nt matches the original PAM nt, don't change it
            #if PAM_TEMPLATE_SEQUENCE[pam_idx] == pam[pam_idx]:
            #    continue
            mutate_nt(pam_idx, pam[pam_idx])

        # save and close all files
        generate_pam_variant_chimera(pam, args.input_pdb, args.output_dir)
