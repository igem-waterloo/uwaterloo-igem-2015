import os

from rosetta import *
from toolbox import mutate_residue
rosetta.init()

def mutate_pose(pose, mutations):
    """Applies list of mutations to the given template pose and returns a mutated version
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
        assert isinstance(aa_num, int)
        assert isinstance(aa_replacement, str) and len(aa_replacement) == 1
        scorefxn = get_fa_scorefxn()
        mutant_pose = mutate_residue(mutant_pose, aa_num, aa_replacement)
    # kims lines from D050 example
    # =================================
    pose_packer = standard_packer_task(mutant_pose)
    pose_packer.restrict_to_repacking()
    # This is a hack, but I want to test. Can't set a movemap, resfiles
    # might be the way to go. Freeze all residues. 
    pose_packer.temporarily_fix_everything()
    # Let's release the PI domain
    for i in range(1110, 1388):
        pose_packer.temporarily_set_pack_residue(i,True)
    packmover = PackRotamersMover(scorefxn, pose_packer)
    packmover.apply(mutant_pose)
    # =================================
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
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)    
    output_pdb_path = os.path.join(output_directory, output_id + ".pdb")
    pose_mutant.dump_pdb(output_pdb_path)
    return output_pdb_path

    
if __name__ == '__main__':
    '''
    This script creates the VQR and EQR mutants from the Kleinstiver 2015 paper
    It is based off the 
    There are difference in residue numbering between Kleinstiver and the
    Anders et al. (2014) experiments that produces the PDB via X-ray
    crystalography. Due to a couple of gaps, and the way PyRosetta handles
    residue numbering, the offset is not constant either. The mutant locations
    in this script have been verified using PyMOL.
    '''
    mutate_pdb("4UN3.tgg.pdb", [(1135+34,'V'),(1335+23,'Q'),(1337+23,'R')],"mutate_folder",'4UN3.VQR')
    mutate_pdb("4UN3.tgg.pdb", [(1135+34,'E'),(1335+23,'Q'),(1337+23,'R')],"mutate_folder",'4UN3.EQR')