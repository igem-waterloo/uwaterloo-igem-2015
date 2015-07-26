import os

import config


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
