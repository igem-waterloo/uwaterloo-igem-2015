from init_genome_camv import init_genome_camv, init_targets_single_P6


# simulation parameters (time in seconds)
complex_concentration = 135000000000
dt = 0.001
t0 = 0.0
t1 = 3600.0
total_turns = int((t1 - t0) / dt)
time_sim = t0

# initialize genome
pseudo_targets = init_targets_single_P6(complex_concentration)
genome_camv = init_genome_camv(pseudo_targets)
genome_camv.initialize_target_cut_probabilities(dt)

# simulation pre-loop behaviour
target_dict = genome_camv.get_targets_from_genome()
open_targets = genome_camv.get_open_targets_from_genome()


for turn in xrange(total_turns):

    # get current targets
    targets_from_genome = genome_camv.get_targets_from_genome()

    # deletion module
    if len(open_targets) > 1:
        # do deletion permutation probabilities for deletion of pairs based on distance

    # cut and repair module
    for key_domain in targets_from_genome.keys():
        domain = genome_camv.domains[key_domain]
        targets_from_domain = domain.targets
        for key_target in targets_from_domain.keys():
            target = targets_from_domain[key_target]
            if target.repaired:  # i.e. not cut
                probability_to_cut = target.cut_probability
                success_cut = ???  # TODO
                if success_cut:
                    target.cut()
                    open_targets.append((key_domain, key_target))

            else:
                probability_to_repair = ???  # TODO
                success_repair = ???  # TODO
                if success_repair:
                    target.repair()
                    open_targets.remove((key_domain, key_target))

    # save turn data (maybe only if stuff happened?)
    # do...

    # update plots if actively showing plots
    # do...

    # increment timer
    time_sim += dt
