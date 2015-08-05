import os

import config
from constants import DNA_ALPHABET


def mutant_variants_dir_by_idx(mutant_idx):
    """Simplify directory retrieval for the pam_variants folder in a specific mutant directory
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
    Returns:
        path to the pam variants directory for mutant_idx
    """
    mutant_variants_dir = os.path.join(config.DIR_PDB_MUTANTS, "mutant_%d" % mutant_idx,
                                       config.FOLDERNAME_MUTANT_VARIANTS)
    return mutant_variants_dir


def mutant_template_dir_by_idx(mutant_idx):
    """Simplify directory retrieval for the template folder in a specific mutant directory
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
    Returns:
        path to the template directory for mutant_idx
    """
    mutant_template_dir = os.path.join(config.DIR_PDB_MUTANTS, "mutant_%d" % mutant_idx,
                                       config.FOLDERNAME_MUTANT_TEMPLATE)
    return mutant_template_dir


def mutant_template_by_idx(mutant_idx):
    """Simplify filename retrieval for a given mutant template pdb based on its index
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
    Returns:
        filename (e.g. "4UN3.mutant_7.pdb")
    """
    return "%s%d.pdb" % (config.MUTANT_TEMPLATE_PREFIX, mutant_idx)


def mutant_template_pdb_path_by_idx(mutant_idx):
    """Simplify mutant pdb path retrieval 
    Args:
        mutant_idx: int >= 0 (0 is default dCas9)
    Returns:
        full filepath to the mutant_idx template pdb, assuming it exists
    Notes:
        - asserts that the path exists because it must point to a file, which shouldn't
          be referenced if it doesn't exist
    """
    mutant_template_pdb_path = os.path.join(mutant_template_dir_by_idx(mutant_idx),
                                            mutant_template_by_idx(mutant_idx))
    assert os.path.exists(mutant_template_pdb_path)
    return mutant_template_pdb_path


def pam_string_from_int(int_pam, length_pam):
    """Given an integer and the length of the pam sequence, return a unique pam
    Args:
        int_pam: integer between 0 and 4 ** (length_pam)
        length_pam: length of the pam sequence (either 3 or 4 currently supported)
    Returns:
        unique 3 or 4 long string of nucleotide characters
    """
    assert 0 <= int_pam < 4 ** length_pam
    if length_pam == 3:
        return DNA_ALPHABET[int_pam / 16] + DNA_ALPHABET[int_pam / 4 % 4] + DNA_ALPHABET[int_pam % 4]
    elif length_pam == 4:
        return DNA_ALPHABET[int_pam / 64] + DNA_ALPHABET[int_pam / 16 % 4] + DNA_ALPHABET[int_pam / 4 % 4] + DNA_ALPHABET[int_pam % 4]
    else:
        raise Exception('Unsupported pam length -- must be 3 or 4 long')
