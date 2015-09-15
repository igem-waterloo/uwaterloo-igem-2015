from init_genome_camv import init_genome_camv, init_targets_all_domains
from random import random
from probabilistic import prob_repair

# simulation parameters (time in seconds)
complex_concentration = 135000000000
dt = 0.1
t0 = 0.0
t1 = 3600.0 # * 18.0
total_turns = int((t1 - t0) / dt)
time_sim = t0

# initialize genome
pseudo_targets = init_targets_all_domains(complex_concentration)
genome_camv = init_genome_camv(pseudo_targets)
genome_camv.initialize_target_cut_probabilities(dt)

# simulation pre-loop behaviour
target_dict = genome_camv.get_targets_from_genome()
open_targets = genome_camv.get_open_targets_from_genome()
probability_to_repair = prob_repair(dt)

# for logging data
data_log = ""
data_file = "simulation_data"


for turn in xrange(total_turns):

    # clear turn log and set to turn
    turn_log = "Time: " + str(turn*dt) + "s\n"

    # get current targets
    targets_from_genome = genome_camv.get_targets_from_genome()

    # deletion module
    #if len(open_targets) > 1:
        # do deletion permutation probabilities for deletion of pairs based on distance
        # regenerate targets_from_genome
        # if large deletion:
        #   first, second = min(target1,target2), max(target1, target2)
        #   turn_log += "Large deletion spanning from", first, "to", second + "\n"

    # cut and repair module
    for key_domain in targets_from_genome.keys():
        domain = genome_camv.domains[key_domain]
        targets_from_domain = domain.targets
        success_cut = False
        success_repair = False
        for key_target in targets_from_domain.keys():
            target = targets_from_domain[key_target]
            if target.repaired:  # i.e. not cut
                probability_to_cut = target.cut_probability
                if random() < probability_to_cut:
                    success_cut = True
                if success_cut:
                    target.cut()
                    open_targets.append((key_domain, key_target))
                    turn_log += target.label + " cut at " + str(target.cut_position) + "\n"

            else:
                if random() < probability_to_repair:
                    success_repair = True
                if success_repair:
                    extra = ""
                    old_sequence = target.sequence
                    old_shift = target.shift
                    target.repair(dt)
                    open_targets.remove((key_domain, key_target))
                    net_indel_size = target.shift - old_shift
                    if old_sequence != target.sequence:
                        extra = "The sequence was changed from " + old_sequence + " to " + target.sequence
                    turn_log += target.label + " repaired at " + str(target.repair_position) + " with an indel of " + str(net_indel_size) + "\n" + extra + "\n"


    # save turn data (maybe only if stuff happened?)
    # \n's count number of events in turn (starts with one)
    if turn_log.count('\n') > 1:
        data_log += turn_log
        print turn_log

    # update plots if actively showing plots
    # do...

    # increment timer
    time_sim += dt

# print data_log
f = open(data_file,'w')
f.write(data_log)
f.close()
