import os

from rosetta import *
from toolbox import mutate_residue


def mutate_pose(pose, mutations):
    """Applies list of mutations to the given pose and returns a pdb
    Args:
        pose: PyRosetta Pose() object representing a loaded pdb structure
        mutations: list of amino acid swaps to apply, format is: [(int, char), ..., (int, char)]
                   where char is a string of length 1 in "ACDEFGHIKLMNPQRSTVWY"
    Returns:
        mutant_pose containing the specified amino acid swaps
    Notes:
        - this procedure doesn't modify the input pose
    """
    mutant_pose = Pose()
    mutant_pose.assign(pose)
    for aa_num, aa_replacement in mutations:
        assert isinstance(aa_num, int)  # could also check that the residue number is an amino acid and in range? note pose residue != cas9 residue
        assert isinstance(aa_replacement, str) and len(aa_replacement) == 1
        mutate_residue(mutant_pose, aa_num, aa_replacement)
    return mutant_pose


def mutate_pdb(input_pdb_path, mutations, output_directory, output_id):
    """Create a new pdb (<output_filename>.pdb) in the output directory containing specified mutations
    Args:
        input_pdb_path: [string] pdb file for template pose (apply mutations to the template)
        mutations: list of amino acid swaps to apply, format is: [(int, char), ..., (int, char)]
                   where str is a character in "ACDEFGHIKLMNPQRSTVWY"
        output_directory: [string] directory to store output pdb file in (e.g. "mutants/some_category/")
        output_id: [string] filename of mutant, do not include ".pdb" (e.g. "some_mutant_id")
    Returns:
        full filepath to the output pdb with the specified mutations
    """
    pose_template = pose_from_pdb(input_pdb_path)
    pose_mutant = mutate_pose(pose_template, mutations)
    if not os.path.exists(output_directory_path):
        os.makedirs(output_directory_path)    
    output_pdb_path = os.path.join(output_directory, output_id + ".pdb")
    mutant_pose.dump(output_pdb_path)
    return output_pdb_path

    
if __name__ == '__main__:
    print "main behaviour not yet implemented"
