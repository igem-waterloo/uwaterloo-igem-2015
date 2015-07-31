import argparse
import errno
import math
import os
import subprocess
import sys

from constants import PAM_TEMPLATE_SEQUENCE, DNA_ALPHABET
from utility import pam_string_from_int


def generate_pam_variant_3dna(pam, template_pdb, output_dir):
    """Mutate base pairs
    Args:
        pam: PAM to mutate to
        template_pdb: input PDB file
        output_dir: output directoty
    """
    mutation = mutation_string(pam)
    template_pdb_basename = os.path.basename(input_pdb)
    new_pdb_path = os.path.join(output_dir, template_pdb_basename[:-4] + "." + pam + ".pdb")
    subprocess.call(["mutate_bases", mutation, input_pdb, new_pdb_path])


def mutation_string(pam):
    """Determine command line output for mutations and positions for call of mutate_base
    Args:
        pam: the PAM to create
    """
    position_idx = { 0 : "5", 1 : "6", 2 : "7", 3 : "8"} # nucleotides we want to mutate are located at positions 5,6,7 in chain D (PAM,NGG) and 8,7,6 in chain C (target strand, NCC).
    position_pairs = { "5" : "8", "6" : "7", "7" : "6", "8" : "5"} # create a dictionary mapping corresponding positions to each other.
    base_pairs = {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'} # create dictionary mapping valid base pairs to each other.
    mutation = ""

    # loop through the PAM sequence and find positions to mutate
    for pam_idx in range(pam_length):
        pos = position_idx[pam_idx]
        base = pam[pam_idx]
        complement_base = base_pairs[base]
        complement_pos = position_pairs[pos]

        # don't bother changing if it matches the input PDB
        if PAM_TEMPLATE_SEQUENCE[pam_idx] == pam[pam_idx]:
            continue
        mutation = mutation + "c=d s=" + pos + " m=D" + base + "; c=c s=" + complement_pos + " m=D" + complement_base + "; "

    mutation = mutation[0:-2] # remove last "; "
    return mutation


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
        generate_pam_variant_3dna(pam_string_from_int(i, pam_length), args.input_pdb, args.output_dir)
