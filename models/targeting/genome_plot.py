import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from math import pi


# old colours
# '#FFE082' sandy
# '#A1887F' light brown
# '#ba68c8' medium pastel purple
# '#779ECB' dark pastel blue
# '#AEC6CF' pastel blue

domain_colours = {'orf': {True: '#81c784',  # pastel green '#80cbc4'
                          False: '#e57373'},  # pastel light red
                  'untracked': {True: '#A1887F',  # light brown
                                False: '#A1887F'},  # light brown
                  'promoter': {True: '#ce93d8',  # pastel light purple
                               False: '#ce93d8'}}  # pastel light purple

target_colours = {'repaired': '#fdfd96',  # pastel yellow
                  'open': 'white',  # just white
                  'inactive': '#fdfd96'}  # pastel yellow

state_colours = {'gene': {'active': '#81c784',
                          'deactivated': '#e57373'},
                 'target': {'cut': '#AEC6CF',
                            'targetable': '#81c784',
                            'untargetable': '#e57373'}}


def genome_plot_polar(genome, genome_label, time=None, output_path=None, flag_show=True):
    # initialize plot
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)

    # plot settings
    ax.grid(False)
    rmax = 20.0
    ax.set_xlim(-rmax, rmax)  # this makes plot look like ring rather than wedges
    ax.set_ylim(-rmax, rmax)  # this makes plot look like ring rather than wedges
    ax.axes.get_xaxis().set_visible(False)  # turn off polar labels
    ax.axes.get_yaxis().set_visible(False)  # turn off polar labels
    theta_init = pi / 2  # where nucleotide zero starts
    theta_direction = -1.0  # clockwise (-1.0), counter-clockwise (+1.0)
    target_radius = -2.0
    domain_radius = 5.0
    text_radius = (domain_radius + 1) * 1.4

    # plot base
    ax.bar(0, domain_radius, width=2*pi, color='lightgray', ec='k', alpha=0.50)

    # prepare data
    genome_length = float(genome.current_length)
    total_targets = 0
    active_targets = 0
    total_genes = 0
    active_genes = 0

    # plot domains and get gene metadata
    target_dict = genome.get_targets_from_genome()
    for domain_key in target_dict.keys():
        # create domain 'patch'
        domain = genome.domains[domain_key]
        color = domain_colours[domain.domain_type][domain.functional]
        if domain.domain_type == 'orf':
            total_genes += 1
            if domain.functional:
                active_genes += 1
        # get domain angular location
        theta_start = (theta_init + theta_direction * 2.0 * pi * domain.domain_start / genome_length) % (2*pi)
        theta_length = (2.0 * pi * (domain.domain_end - domain.domain_start) / genome_length) % (2*pi)
        theta_mid = theta_start + theta_direction * theta_length / 2
        ax.bar(theta_start, domain_radius, width=theta_direction*theta_length, color=color, ec='k', alpha=0.75)
        ax.plot(theta_mid, domain_radius, 'o', color=color)
        #ax.plot(theta_mid, 0, 'o', color=color)
        # write domain label
        text_rotation_degrees = ((theta_mid * 180 / pi) + 270) % 360
        plt.text(theta_mid, text_radius, "%s:\n%s" % (domain.domain_type, domain.label), fontsize=10,
                 rotation=text_rotation_degrees, horizontalalignment='center', verticalalignment='center')

    # plot targets and get target metadata
    for domain_key in target_dict.keys():
        domain = genome.domains[domain_key]
        targets_from_domain = domain.targets
        for target_key in targets_from_domain.keys():
            total_targets += 1
            target = targets_from_domain[target_key]
            # get target angular location
            phi_direction = theta_direction * target.sense
            phi_start = (theta_init + theta_direction * 2.0 * pi * target.current_start / genome_length) % (2*pi)
            phi_length = (2.0 * pi * len(target.sequence) / genome_length) % (2*pi)
            phi_mid = phi_start + phi_direction * phi_length / 2
            target_rotation_degrees = ((phi_mid * 180 / pi) + 270) % 360
            plt.text(phi_mid, target_radius, "Cuts: %d" % target.total_cuts, fontsize=10,
                     rotation=target_rotation_degrees, horizontalalignment='center', verticalalignment='center')
            if target.repaired:
                color_target = target_colours['repaired']
                ax.bar(phi_start, domain_radius, width=phi_direction*phi_length, color=color_target, ec='k', alpha=0.75)
            else:
                color_target = target_colours['open']
                open_buffer = 1.0
                eta_start = phi_start - phi_direction*phi_length*open_buffer
                # plot big white bar
                ax.bar(eta_start, domain_radius, width=phi_direction*phi_length*(1+2*open_buffer),
                       color=color_target, ec='white', alpha=1.0)
                # add edges back
                ax.bar(eta_start, domain_radius, width=0, color=color_target, ec='k', alpha=0.75)
                ax.bar(eta_start + phi_direction*phi_length*(1+2*open_buffer), domain_radius, width=0,
                       color=color_target, ec='k', alpha=0.75)
            if target.targetable:  # TODO FIND WAY TO ANNOTATE
                active_targets += 1
                #ax.plot(phi_mid, target_radius, '+', color='k', markersize=12, markeredgewidth=2)
            else:
                #ax.plot(phi_mid, target_radius, '_', color='k', markersize=12, markeredgewidth=2)
                continue

    # plot genome metadata
    fs_title = 18
    fs_data = 14
    ax.annotate("Genome: %s" % genome_label, xy=(2, 2), xytext=(pi/2, rmax * 0.65), textcoords='data', fontsize=fs_title,
                horizontalalignment='center', verticalalignment='center')
    ax.annotate("Length: %d" % int(genome_length), xy=(2, 2), xytext=(pi/2, -rmax*0.8), textcoords='data', fontsize=fs_data,
                horizontalalignment='center', verticalalignment='center')
    ax.annotate("Active Targets: %d/%d" % (active_targets, total_targets), xy=(2, 2), xytext=(pi/2, -rmax*0.95), textcoords='data', fontsize=fs_data,
                horizontalalignment='center', verticalalignment='center')
    ax.annotate("Functional Genes: %d/%d" % (active_genes, total_genes), xy=(2, 2), xytext=(3*pi/2, -rmax*0.9), textcoords='data', fontsize=fs_data,
                horizontalalignment='center', verticalalignment='center')
    if time is not None:
        ax.annotate("Time: %.2f min" % time, xy=(2, 2), xytext=(-pi/2, rmax*0.75), textcoords='data', fontsize=fs_data,
                    horizontalalignment='center', verticalalignment='center')

    if output_path is not None:
        fig.set_size_inches(20.0, 8.0)  # alternative: 20.0, 8.0
        fig.tight_layout()
        plt.savefig(output_path)

    if flag_show:
        plt.show()

    return fig


def plot_states(states, datatype, labels_to_plot=None, output_path=None, flag_show=True):
    """Plots timeseries of statedata as horizontal coloured bars
    Args:
        states: dictionary of label: state vector
        datatype: either "gene" or "target"
        labels_to_plot: [default:None] ordered list that's a subset of the keys from states; plot all if None
    """
    assert datatype in ['gene', 'target']
    colour_dict = state_colours[datatype]
    data_labels = states.keys()
    time = [float(t) for t in states['time']]
    if labels_to_plot is None:
        labels_to_plot = data_labels
        labels_to_plot.remove('time')
        labels_to_plot.remove('header')
        labels_to_plot.sort()
    length_data = len(time)
    length_labels = len(labels_to_plot)

    fig = plt.figure()
    ax = plt.gca()
    x0 = time[0]
    y0 = 0.0
    dx = time[1] - time[0]
    dy = dx
    x = x0
    y = y0
    for i, label in enumerate(labels_to_plot):
        y = y0 + dy*i
        state_data = states[label]
        for j, elem in enumerate(state_data):
            x = x0 + dx*j
            ax.add_patch(mpatches.Rectangle((x - dx*0.5, y), width=dx, height=dy*0.5, color=colour_dict[elem])) #, ec='k'))
    ax.set_xlabel('time')
    ax.set_ylabel(datatype)
    ax.set_xticks([time[i] for i in xrange(0, len(time), len(time) / 20)])  # downsample x-axis ticks
    ax.set_yticks([y0 + (i+0.25)*dy for i in xrange(length_labels)])
    ax.set_yticklabels(labels_to_plot)
    ax.set_xlim(x0 - 0.5*dx, (length_data - 1)*dx*1.05)
    ax.set_ylim(y0 - 0.3*dy, (length_labels)*dy)

    if output_path is not None:
        fig.set_size_inches(20.0, 8.0)  # alternative: 20.0, 8.0
        fig.tight_layout()
        plt.savefig(output_path)

    if flag_show:
        plt.show()

    return


def plot_state_multirun(states, datatype, labels_to_plot=None, output_path=None, flag_show=True):
    """Plots timeseries of gene or target state (average) condition
    Args:
        states: dictionary of label: state vector
        datatype: either "gene" or "target"
        labels_to_plot: [default:None] ordered list that's a subset of the keys from states; plot all if None
    """
    assert datatype in ['gene', 'target']
    data_labels = states.keys()
    time = [float(t) for t in states['time']]
    if labels_to_plot is None:
        labels_to_plot = data_labels
        labels_to_plot.remove('time')
        labels_to_plot.remove('header')
        labels_to_plot.sort()
    length_data = len(time)
    length_labels = len(labels_to_plot)
    # plot
    fig = plt.figure()
    ax = plt.gca()
    for i, label in enumerate(labels_to_plot):
        state_data = states[label]
        plt.plot(time, state_data)
    ax.set_xlabel('time')
    ax.set_ylabel(datatype + ' fraction active')
    """
    ax.set_xticks([time[i] for i in xrange(0, len(time), len(time) / 20)])  # downsample x-axis ticks
    ax.set_xlim(x0 - 0.5*dx, (length_data - 1)*dx*1.05)
    """
    # output plot
    if output_path is not None:
        fig.set_size_inches(20.0, 8.0)  # alternative: 20.0, 8.0
        fig.tight_layout()
        plt.savefig(output_path)
    if flag_show:
        plt.show()
    return


if __name__ == '__main__':
    from init_genome_camv import init_genome_camv, init_targets_all_domains
    complex_concentration = 22.101  # nM
    dt = 0.1
    pseudo_targets = init_targets_all_domains(complex_concentration)
    genome_camv = init_genome_camv(pseudo_targets)
    genome_camv.initialize_target_cut_probabilities(dt)
    print "Generating test plot..."
    genome_plot_polar(genome_camv, 'CaMV', time=60.0, output_path='test_genome.png', flag_show=True)
    print "Done"
