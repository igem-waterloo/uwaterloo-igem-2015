from rosetta import *
import os

gap_sequence = 'YEKLKGSPEDN'
'''
We need to add amino acids to the end and beginning because PyRosetta will not
allow N-terminal or C-terminal AAs to be appended into the sequences because
REASONS
'''
graft_pose_sequence = 'A' + gap_sequence + 'W'
graft_pose = Pose()
graft_pose = pose_from_sequence(graft_pose_sequence, 'fa_standard', True)

cas9_gapped = pose_from_pdb('4UN3.tgg.pdb')
pose_gap_start = 1276
insert_length = len(gap_sequence)

for i in range(0,insert_length):
    # We start from index 0, want the 2nd graft_pose residue
    residue = graft_pose.residue(i+2)
    # pose_gap_start is first 'gap' residue, seqpos is the residue to 
    # append the new one to, so decrease by one
    seqpos = pose_gap_start + i - 1
    cas9_gapped.append_polymer_residue_after_seqpos(residue, seqpos, False)

# Create a Fragment mover to get the residues into place. Does not work, maybe
# not proper fragment file?
fragset = ConstantLengthFragSet(3)
fragfile_path = os.path.abspath("C:\\Program Files\\PyRosetta\\test\\data\\workshops\\aat000_03_05.200_v1_3")
fragset.read_fragment_file(fragfile_path)
movemap = MoveMap()
movemap.set_bb_true_range(pose_gap_start,pose_gap_start+insert_length-1)
mover_3mer = ClassicFragmentMover(fragset, movemap)

# Try small and shear movers. Limits to the proper residues, but will move all
# over the place. Try with Monte Carlo acceptance?
kT = 1.0
n_moves = 100
small_mover = SmallMover(movemap, kT, n_moves)
shear_mover = ShearMover(movemap, kT, n_moves)
small_mover.apply(cas9_gapped)
shear_mover.apply(cas9_gapped)

# Try Minimization
minmover = MinMover()
minmover.movemap(movemap)
minmover.score_function(get_fa_scorefxn())
