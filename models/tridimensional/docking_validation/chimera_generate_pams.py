import os
import argparse
import chimera
from chimera import runCommand
from chimera import replyobj

dna_nts = ['a','c','g','t']


def mutate_nt(pam_idx, base):
    """Mutate a base pair at pam_idx to base
    Args:
        pam_idx: the index of the PAM nucleotide to be modified
        base: what the nt at pam_idx will be changed to
    """
    position_idx = { 0 : "5.d", 1 : "6.d", 2 : "7.d"} # Nucleotides we want to mutate are located at positions 5,6,7 in chain D (PAM,NGG) and 8,7,6 in chain C (target strand, NCC).
    position_pairs = { "5.d" : "8.c", "6.d" : "7.c", "7.d" : "6.c"} # Create a dictionary mapping corresponding positions to each other.
    base_pairs = {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'} # Create dictionary mapping valid base pairs to each other.

    pos = position_idx[pam_idx]
    complement_base = base_pairs[base]
    complement_pos = position_pairs[pos]
    runCommand("swapna " + base + " : " + pos )
    runCommand("swapna " + complement_base + " : " + complement_pos)

def generate_pdb(original_pdb_basename, pam_seq):
    """Generate a new PDB file
    Args:
        original_pdb_basename: original fine basename
        pam_seq: new PAM sequence to be listed in name
    """
    new_pdb_path = os.path.join(args.output_dir, original_pdb_basename[:4] + "." + pam_seq + ".pdb")
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

    assert args.num_pams
    assert args.input_pdb
    assert args.output_dir
    args.num_pams = int(args.num_pams) # convert string from command line to int
    assert 64 == args.num_pams or 256 == args.num_pams
    assert os.path.isfile(args.input_pdb)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

for i in range(args.num_pams):
    if 64 == args.num_pams:
        pam = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]
    else:
        pam = dna_nts[i / 64] + dna_nts[i / 16 % 4] + dna_nts[i / 4 % 4] + dna_nts[i % 4]

    # TGG PAM site already exists, let's skip it
    if "tgg" == pam:
        continue
    # open up the file again each time, for now
    runCommand("open " + args.input_pdb)
    # Loop through the PAM sequence and mutate positions
    for n in range(len(pam)):
        base = pam[n]
        # If the first base is T, don't mutate
        if base == 't' and n == 0:
            continue
        # If the 2nd or third bases are G, don't mutate
        if all([base == 'g', any([n == 1, n == 2])]):
            continue
        mutate_nt(n, base)
    # Save and close all files
    generate_pdb(os.path.basename(args.input_pdb), pam)
