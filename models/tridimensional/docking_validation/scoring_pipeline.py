import os

from mutate_pdb import mutate_pdb
from results_csv import results_to_csv
from config import MUTANT_DATABASE_PATH, MUTANT_TEMPLATE_PREFIX, PATH_PDB_SPECIAL_ORIGINAL
from utility import mutant_variants_dir_by_idx, mutant_template_dir_by_idx, mutant_template_by_idx, mutant_template_pdb_path_by_idx, 


# TODO
"""
A
mutation_list = ...

B
pam_variants_chimera(num_pam_variants, mutant_template_pdb_path, mutant_variants_pdb_directory)
# or
pam_variants_3dna(num_pam_variants, mutant_template_pdb_path, mutant_variants_pdb_directory)
"""


def get_mutations_from_index(mutant_idx):
    """Given an index, extract the corresponding mutations from the database
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
        mutant_database_path: location of mutant information file
    Returns:
        list of mutations of the form [(aa_num, aa_char),...]
    """
    mutation_list = ...
    return mutation_list


def generate_mutant_from_index(mutant_idx, template_pdb_path=SPECIAL_PDB_ORIGINAL):
    """Given an index, make the mutant with the corresponding mutation code
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
        template_pdb_path: full path (folder and filename) of the template pdb
    Returns:
        path to the mutant pdb
    """
    mutation_list = get_mutations_from_index(mutant_idx)
    mutant_template_dir = mutant_template_dir_by_idx(mutant_idx)
    mutant_template_id = "%s%d" % (MUTANT_TEMPLATE_PREFIX, mutant_idx)
    mutant_pdb_path = mutate_pdb(template_pdb_path, mutation_list, mutant_template_dir, mutant_template_id)
    return mutant_pdb_path


def generate_pam_variants_from_mutant(mutant_idx, num_pam_variants, mutant_pdb_subdirectory):
    """Generate all the pam variants for a given mutant 
    Args:
        mutant_idx:  int >= 0 (0 is default dCas9)
        num_pam_variants: 64 or 256 -- represents all 3 or 4nt pam variants
        mutant_pdb_directory: folder of name "mutant_INT" where int >= 0
    Returns:
        path to the pam variant folder (containing 64 or 256 pdbs)
    """
    assert num_pams == 64 or num_pams == 256  # currently support 3 or 4 nt PAM sites
    mutant_template_pdb_path = mutant_template_pdb_path_by_idx(mutant_idx)
    mutant_variants_pdb_dir = mutant_variants_dir_by_idx(mutant_idx) 
    # ========================================================================
    pam_variants_chimera(num_pam_variants, mutant_template_pdb_path, mutant_variants_pdb_dir)
    # or
    pam_variants_3dna(num_pam_variants, mutant_template_pdb_path, mutant_variants_pdb_dir)
    # ========================================================================
    return mutant_variants_pdb_dir


if __name__ == '__main__':
    print "Main not implemented"
    # get mutant int from args
    # get k (dock 256 run repeats) from args
