import argparse
import datetime
import os
from time import time

from rosetta import *

from results_csv import results_to_csv


# variations to account for
dna_nts = "acgt"
programs = ["Chimera", "3DNA"]


def write_dock_stats(score_directory, filename, dock_stats, time_diff):
    """Writes separate files for each variant
    Notes:
    - score_directory has the form results/timestamped folder/
    - dock stats has the form: [dna_init, dna_final, fa_init, fa_final]
    - 'a' appends in case of accidental duplication
    - scores and times are passed in explicitly for argument clarity
    """
    path = os.path.join(score_directory, filename)
    f = open(path, 'a')
    f.write("Initial DNA score: %8.3f\n" % dock_stats[0])
    f.write("Final DNA Score: %8.3f\n" % dock_stats[1])
    f.write("Initial FA score: %8.3f\n" % dock_stats[2])
    f.write("Final FA score: %8.3f\n" % dock_stats[3])
    f.write("Total variant time: %8.3f\n" % time_diff)
    f.close()
    return


def dock_simple(pose):
    """Coarse docking of a pose representing a PAM / program variant
    Returns:
        list of scores in the form [dna_init, dna_final, fa_init, fa_final]
    Notes:
    - docking optimizes for DNA score, which is weighted differently than "Full Atom" (FA) score
    - dna score from demo/D110_DNA_interface.py
    Potential Bugs:
    - setup_foldtree(...) crashes on some systems/configurations
    - docking.set_partners(...) may not be needed
    """
    # specify foldtree for simple docking
    setup_foldtree(pose, 'B_CD', Vector1([1]))

    # specify scoring functions
    fa_score = get_fa_scorefxn()  # standard full atom score
    dna_score = create_score_function('dna')  # 
    dna_score.set_weight(fa_elec, 1)

    # specify docking protocol
    docking = DockMCMProtocol()
    docking.set_scorefxn(dna_score)
    docking.set_scorefxn_pack(fa_score)
    docking.set_partners("B_ACD")

    # obtain initial and final scores after docking
    dna_init = dna_score(pose)
    fa_init = fa_score(pose)
    docking.apply(pose)
    dna_final = dna_score(pose)
    fa_final = fa_score(pose)
    
    return [dna_init, dna_final, fa_init, fa_final]


def dock_variants(pam_variants, path_to_variant_scores):
    """Docks and scores 2 pdbs for each PAM variant (one for each nt program) using simple docking
    Args:
        pam_variants: list of integers (any from 0 to 63 without repeats) which map to pam strings
        path_to_variant_scores: path to the subdirectory of "results" where the variants are stored
    Notes:
    - creates a text file (e.g. 'results_agg_Chimera.txt') for each variant
    - path to variants is typically root/results/<timestamped folder>/<variants>
    - assumes current directory is the root of a folder that contains pdbs in Chimera and 3DNA directories
    """
    for i in pam_variants:
        variant = dna_nts[i / 16] + dna_nts[i / 4 % 4] + dna_nts[i % 4]
        for program in programs:
            # track runtime while loading and passing pose to the simple docker
            print "Running for variant: %s_%s" % (variant, program)
            time_init = time()
            pdb_filename = "4UN3." + variant + ".pdb"
            loaded_pose = pose_from_pdb(program + os.sep + pdb_filename)
            dock_stats = dock_simple(loaded_pose)
            time_final = time()

            # write results to file
            results_filename = variant + "_" + program + ".txt"
            write_dock_stats(path_to_variant_scores, results_filename, dock_stats, time_final - time_init)
            print "Finished writing scores for variant: %s_%s" % (variant, program)
    return


if __name__ == '__main__':
    # parse arguments 
    parser = argparse.ArgumentParser(description='Run scoring on PAM variants')
    parser.add_argument('-s', '--start', metavar='N', type=int, help='Starting PAM number (0 = aaa)')
    parser.add_argument('-e', '--end', metavar='N', type=int, help='Ending PAM number (63 = ttt)')
    args = parser.parse_args()
    assert 0 <= args.start < args.end <= 63
    pam_variants_to_score = range(args.start, args.end + 1)

    # specify location of pdb directories
    #os.chdir("")  # specify pdb folder root if script isn't located there

    # setup results directory and timestamped variant score subdirectory
    results_folder = "results"
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
    path_to_variant_scores = "results" + os.sep + timestamp + os.sep
    os.makedirs(path_to_variant_scores)
    
    # initialize pyrosetta and score variants
    init(extra_options="-mute all")  # reduce rosetta print calls
    dock_variants(pam_variants_to_score, path_variants)

    # collect text files into a csv
    results_to_csv(path_to_variant_scores)
