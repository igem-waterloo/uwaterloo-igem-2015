from rosetta import *

init()

pose = pose_from_pdb("putsomefilehere")

pyobs = PyMOL_Observer()
pyobs.add_observer(pose)
pyobs.pymol.keep_history(True)
# pyobs.generalEvent()      // ???

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