import os
from rosetta import *
init()

# Define 'center of particles' given list of 3 tuples with xyz coordinates
def findParticlesCenter(residue,list_of_atoms):
    xyz = [0,0,0]
    length = len(list_of_atoms)
    for atom in list_of_atoms:
        xyz[0] += residue.xyz(atom)[0]
        xyz[1] += residue.xyz(atom)[1]
        xyz[2] += residue.xyz(atom)[2]
    xyz[:] = [x / length for x in xyz]
    return xyz

def vectorDiff(xyz_1,xyz_2):
    xyz_shift = [x_1 - x_2 for x_1, x_2 in zip(xyz_1,xyz_2)]
    return xyz_shift

# Return a center of particles for backbone, main chain and side chain for a given residue
# This could probably be abstracted to take in a residue and a list of lists of atoms, and 
# return positions as a list
def findMajor3Pos(residue):
    backbone = range(1,residue.last_backbone_atom()+1)
    mainchain = residue.mainchain_atoms()
    sidechain = range(residue.first_sidechain_atom(),residue.attached_H_begin(1))

    backbone_pos = findParticlesCenter(residue,backbone)
    mainchain_pos = findParticlesCenter(residue,mainchain)
    sidechain_pos = findParticlesCenter(residue,sidechain)

    major3Pos = [backbone_pos, mainchain_pos, sidechain_pos]
    return major3Pos

# Find changes in atom positions
def findResidueDiffs(residue1,residue2):
    residue1_pos = findMajor3Pos(residue1)
    residue2_pos = findMajor3Pos(residue2)
    shifts = []
    for j in range(len(residue1_pos)):
        shifts.append(vectorDiff(residue1_pos[j],residue2_pos[j]))
        shifts[j][:] = [ round(elem, 5) for elem in shifts[j]]
    return shifts

dna_nts = "acgt"
programs = ["Chimera", "3DNA"]
variants =[]
residue_list = [1393,1394,1395,1419,1420,1421]
resultFile = "diffresults.txt"

for k in range(0, 64):
    # Creates a list of three nucleotide PAM sequences.
    variants.append(dna_nts[k / 16] + dna_nts[k / 4 % 4] + dna_nts[k % 4])

for pam in variants:
    chimeraPose = pose_from_pdb("Chimera/4UN3." + pam + ".pdb")
    dna3Pose = pose_from_pdb("3DNA/4UN3." + pam + ".pdb")

    for resi in residue_list:
        resi3DNA = dna3Pose.residue(resi)
        resiChimera = chimeraPose.residue(resi)
        shifts = findResidueDiffs(resiChimera,resi3DNA)
        with open(resultFile, "a") as f:
            f.write(pam + ",")
            f.write(str(resi) + ",")
            for i in range(len(shifts)):
                for k in range(len(shifts[i])):
                    f.write(str(shifts[i][k]) + ",")
            f.write("\n")
            f.close()
