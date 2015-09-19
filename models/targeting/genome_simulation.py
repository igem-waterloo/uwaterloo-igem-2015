import datetime
import matplotlib.pyplot as plt
import os
import random

import make_video
from genome_csv import results_to_csv, csv_to_dict, map_genome_events, map_target_events
from genome_plot import genome_plot_polar, plot_states
from init_genome_camv import init_genome_camv, init_targets_all_domains, init_targets_multi_P6
from probabilistic import prob_repair


def genome_simulate(flag_plot=True, flag_multirun=False, batch_data_path=None):

    # output management
    if not flag_multirun:
        runs_folder = "runs" + os.sep  # store timestamped runs here
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
        time_folder = current_time + os.sep
        current_run_folder = runs_folder + time_folder
        # subfolders in the timestamped run directory:
        data_folder = os.path.join(current_run_folder, "data")
        plot_genome_folder = os.path.join(current_run_folder, "plot_genome")
        plot_data_folder = os.path.join(current_run_folder, "plot_data")
        # create dirs conditionally
        dir_list = [runs_folder, current_run_folder, data_folder]
        if flag_plot:
            dir_list += [plot_genome_folder, plot_data_folder]
        for dirs in dir_list:
            if not os.path.exists(dirs):
                os.makedirs(dirs)
    else:
        assert batch_data_path is not None
        assert not flag_plot  # unfortunately, we don't support plotting during batch runs
        data_folder = batch_data_path
        # create dirs conditionally
        dir_list = [data_folder]
        for dirs in dir_list:
            if not os.path.exists(dirs):
                os.makedirs(dirs)

    # simulation parameters (time in seconds)
    complex_concentration = 22.101  # nM
    dt = 1.0
    t0 = 0.0
    t1 = 2.0 * 18.0 * 3600.0
    total_turns = int((t1 - t0) / dt)
    time_sim = t0
    plot_period = 30  # in turns
    plot_count = 0

    # initialize genome
    pseudo_targets = init_targets_multi_P6(complex_concentration)
    genome_camv = init_genome_camv(pseudo_targets)
    genome_camv.initialize_target_cut_probabilities(dt)

    # simulation pre-loop behaviour
    target_dict = genome_camv.get_targets_from_genome()
    open_targets = genome_camv.get_open_targets_from_genome()
    probability_to_repair = prob_repair(dt)
    double_cut_probability = 1.55*10.0**(-5)  # see Tessa

    # for logging data
    data_log = ""
    data_file = os.path.join(data_folder, "simulation_data.txt")

    # variables for csv writing
    genome_events = [0]*total_turns
    target_events = [0]*total_turns
    genome_header = ["time"] + genome_camv.domains.keys()
    target_header = ["time"]
    for key in genome_camv.domains.keys():
        target_header += genome_camv.get_targets_from_genome()[key].keys()

    def target_state(target_dict):
        if {} == target_dict:
            return None
        for key, value in target_dict.iteritems():
            if not value.repaired:
                return "cut"
            elif value.targetable:
                return "targetable"
            else:
                return "untargetable"

    for turn in xrange(total_turns):

        # clear turn log and set to turn
        turn_log = "Time: " + str(turn*dt) + "s\n"

        # get current targets
        targets_from_genome = genome_camv.get_targets_from_genome()
        open_targets = genome_camv.get_open_targets_from_genome()

        # place data in rows for csv to write later
        genome_events[turn] = map_genome_events(str(turn*dt), genome_camv.domains, genome_header)
        target_events[turn] = map_target_events(str(turn*dt), targets_from_genome, target_header)

        # deletion module
        if len(open_targets) > 1:
            # time_with_double_cut += dt
            double_cut_success = False
            if random.random() < double_cut_probability:
                double_cut_success = True
            if double_cut_success:
                targets = random.sample(open_targets, 2)
                target1 = targets_from_genome[targets[0][0]][targets[0][1]]
                target2 = targets_from_genome[targets[1][0]][targets[1][1]]
                first = min(target1.current_start, target2.current_start)
                second = max(target1.current_start, target2.current_start)
                genome_camv.large_deletion(target1, target2, dt)
                turn_log += "Large deletion spanning from " + str(first) + " to " + str(second) + "\n"
                targets_from_genome = genome_camv.get_targets_from_genome()

        # cut and repair module
        for key_domain in targets_from_genome.keys():
            domain = genome_camv.domains[key_domain]
            targets_from_domain = domain.targets
            for key_target in targets_from_domain.keys():
                success_cut = False
                success_repair = False
                target = targets_from_domain[key_target]
                if target.repaired:  # i.e. not cut
                    probability_to_cut = target.cut_probability
                    if random.random() < probability_to_cut:
                        success_cut = True
                    if success_cut:
                        target.cut()
                        open_targets.append((key_domain, key_target))
                        turn_log += target.label + " cut at " + str(target.cut_position) + "\n"

                else:
                    if random.random() < probability_to_repair:
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
            if not flag_multirun:
                print turn_log

        # update plots if actively showing plots
        if turn % plot_period == 0 and flag_plot:
            plot_path = os.path.join(plot_genome_folder, "genome_%05d.png" % plot_count)
            genome_plot_polar(genome_camv, 'CaMV', time=time_sim/60.0, output_path=plot_path, flag_show=False)
            plt.close()
            plot_count += 1

        # increment timer
        time_sim += dt

    # write data to csvs
    csv_states_gene = "states_gene.csv"
    csv_states_target = "states_target.csv"
    results_to_csv(data_folder, csv_states_gene, genome_header, genome_events)
    results_to_csv(data_folder, csv_states_target, target_header, target_events)

    # print data_log
    f = open(data_file, 'w')
    f.write(data_log)
    f.close()

    # create data plots
    if flag_plot:
        states_gene = csv_to_dict(os.path.join(data_folder, csv_states_gene))
        states_target = csv_to_dict(os.path.join(data_folder, csv_states_target))
        domains_to_plot = list(set([elem["domain_label"] for elem in pseudo_targets]))
        plot_states(states_gene, "gene", labels_to_plot=domains_to_plot, output_path=os.path.join(plot_data_folder, "states_gene.png"), flag_show=False)
        plot_states(states_target, "target", output_path=os.path.join(plot_data_folder, "states_target.png"), flag_show=False)

    # create video of genome plots
    if flag_plot:
        fps = 15
        video_path = os.path.join(current_run_folder, "genome_%dmin_%dfps.mp4" % (int(t1/60.0), fps))
        make_video.make_video_ffmpeg(plot_genome_folder, video_path, fps=fps)

    return


if __name__ == '__main__':
    genome_simulate(False)
