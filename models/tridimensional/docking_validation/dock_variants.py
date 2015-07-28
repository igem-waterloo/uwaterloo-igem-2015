import argparse
import datetime
import errno
import os
from time import time

from rosetta import *

from constants import SCOREFILE_LINES
from results_csv import results_to_csv
from utility import pam_string_from_int


# variations to account for
programs = ["Chimera", "3DNA"]


def write_dock_stats(score_directory, filename, dock_stats, time_diff_total, time_diff_docking):
    """Writes separate files for each variant
    Notes:
    - score_directory has the form results/timestamped folder/
    - dock stats has the form: [fa_init, fa_final, dna_init, dna_final]
    - 'a' appends in case of accidental duplication
    - scores and times are passed in explicitly for argument clarity
    """
    path = os.path.join(score_directory, filename)
    f = open(path, 'a')
    f.write("Initial FA score: %8.3f\n" % dock_stats[0])
    f.write("Final FA score: %8.3f\n" % dock_stats[1])
    f.write("Initial DNA score: %8.3f\n" % dock_stats[2])
    f.write("Final DNA Score: %8.3f\n" % dock_stats[3])
    f.write("Total variant time: %8.3f\n" % time_diff_total)
    f.write("Docking variant time: %8.3f\n" % time_diff_docking)
    f.close()
    return


def dock_simple(pose):
    """Coarse docking of a pose representing a PAM / program variant
    Returns:
        list of scores in the form [fa_init, fa_final, dna_init, dna_final]
    Notes:
    - docking optimizes for DNA score, which is weighted differently than "Full Atom" (FA) score
    """
    # specify scoring functions
    fa_score = get_fa_scorefxn()
    dna_score = create_score_function('dna')
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
    return [fa_init, fa_final, dna_init, dna_final]


def dock_complex(pose):
    """Complex docking of a pose representing a PAM / program variant
    Returns:
        list of scores in the form [fa_init, fa_final, dna_init, dna_final]
    Notes:
    - not implemented
    """
    raise Exception("Complex docking not implemented")


def dock_variants(pam_variants, path_to_scores, path_to_pdbs='', complex_docking_flag=False, pam_length=4):
    """Docks and scores 2 pdbs for each PAM variant (one for each nt program) using simple docking
    Args:
        pam_variants: list of integers (any from 0 to 63 without repeats) which map to pam strings
        path_to_scores: path to the subdirectory of "results" where the variants are stored
        path_to_pdbs: [default: current directory] path to location of chimera/3DNA folders of PAM variants
        complex_docking_flag: [default: False] if True, use complex dock function (NOT IMPLEMENTED)
        pam_length: [default: 4] length of the pam sequence to be investigated
    Notes:
    - creates a text file (e.g. 'results_agg_Chimera.txt') for each variant
    - path to variants is typically root/results/<timestamped folder>/<variants>
    - assumes current directory is the root of a folder that contains pdbs in Chimera and 3DNA directories
    """
    for idx in pam_variants:
        variant = pam_string_from_int(idx, pam_length)
        for program in programs:
            print "Running for variant: %s_%s" % (variant, program)
            pdb_path = os.path.join(path_to_pdbs, program, "4UN3." + variant + ".pdb")

            # track runtime while loading and passing pose to the simple docker
            time_init_total = time()
            loaded_pose = pose_from_pdb(pdb_path)

            time_init_docking = time()
            if complex_docking_flag:
                dock_stats = dock_complex(loaded_pose)
            else:
                dock_stats = dock_simple(loaded_pose)

            time_final = time()
            time_diff_total = time_final - time_init_total
            time_diff_docking = time_final - time_init_docking

            # write results to file
            results_filename = variant + "_" + program + ".txt"
            write_dock_stats(path_to_scores, results_filename, dock_stats, time_diff_total, time_diff_docking)
            print "Finished writing scores for variant: %s_%s" % (variant, program)
    return


if __name__ == '__main__':
    # create parser and parse arguments
    parser = argparse.ArgumentParser(description='Run docking and scoring on PAM variants')
    parser.add_argument('-s', '--start', metavar='N', type=int,
                        help='starting PAM number (0 = aaaa), inclusive')
    parser.add_argument('-e', '--end', metavar='N', type=int,
                        help='ending PAM number (255 = tttt), inclusive')
    parser.add_argument('-o', '--output_dir', metavar='D', type=str,
                        help='path to output directory for variant scores')
    parser.add_argument('--pdb_dir', metavar='D', default='', type=str,
                        help='path to root of variant pdb directories (default: "")')
    parser.add_argument('--complex', metavar='B', nargs='?', const=True, default=False,
                        type=str, help='[switch] select complex docking (default: "False")')
    parser.add_argument('--csv', metavar='B', nargs='?', const=True, default=False,
                        type=str, help='[switch] compile scores to csv (default: "False")')
    parser.add_argument('--pam64', metavar='B', nargs='?', const=3, default=4,
                        type=str, help='[switch] assume 64 pam variants (default: assume 256)')
    args = parser.parse_args()

    # setup range of pam variants
    assert 0 <= args.start <= args.end <= (4 ** args.pam64) - 1
    pam_variants_to_score = range(args.start, args.end + 1)

    # setup output path for scoring
    if args.output_dir is not None:
        path_to_scores = args.output_dir
    else:
        results_folder = "results"
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
        path_to_scores = os.path.join("results", timestamp)
    try: # check existence again to handle concurrency problems
        os.makedirs(path_to_scores)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path_to_scores):
            pass
        else:
            raise

    # initialize pyrosetta and score variants
    init(extra_options="-mute all")  # reduce rosetta print calls
    dock_variants(pam_variants_to_score, path_to_scores,
                  path_to_pdbs=args.pdb_dir, complex_docking_flag=args.complex, pam_length=args.pam64)

    # collect score txt files into a csv
    if args.csv:
        results_to_csv(path_to_scores)
