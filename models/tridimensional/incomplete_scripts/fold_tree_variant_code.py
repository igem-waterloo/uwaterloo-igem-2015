# Manual Fold Tree Setup
        original_pose = pose_from_pdb(program + "/4UN3." + variant + ".pdb")    
        pose_fold_tree = FoldTree(original_pose.total_residue())
        pose_fold_tree.new_jump(1, 1388, 1387)
        pose_original.fold_tree(pose_fold_tree)

# Default Fold Tree
# Blank

# BC Fold Tree Setup
        setup_foldtree(original_pose, 'B_C', Vector1([1]))

# BD Fold Tree Setup
        setup_foldtree(original_pose, 'B_D', Vector1([1]))