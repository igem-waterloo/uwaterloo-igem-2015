import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from math import pi


domain_colours = {'orf': '#66FF66',
                  'ncr': '#9999FF',
                  'promoter': '#FF9999'}


def genome_plot_polar(genome):
    # initialize plot
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    patches = []

    # plot settings
    ax.grid(False)
    #ax.set_rmax(20.25)   # set maximum radius
    rmax = 20.0
    ax.set_xlim(-rmax, rmax)  # this makes plot look like ring rather than wedges
    ax.set_ylim(-rmax, rmax)  # this makes plot look like ring rather than wedges
    ax.axes.get_xaxis().set_visible(False)  # turn off polar labels
    ax.axes.get_yaxis().set_visible(False)  # turn off polar labels
    theta_init = pi / 2  # where does nucleotide zero start
    theta_direction = -1.0  # clockwise (-1.0), counter-clockwise (+1.0)
    domain_radius = 5.0
    domain_width = 1.0
    text_radius = (domain_radius + 1) * 1.4

    # plot genome_metadata
    genome_length = float(genome.length)
    ax.bar(0, domain_radius, width=2*pi, color='lightgray', ec='k', alpha=0.50)
    ax.annotate("Genome: CaMV", xy=(2,2), xytext=(pi/2, rmax * 0.65), textcoords='data', fontsize=16, horizontalalignment='center', verticalalignment='center')
    ax.annotate("Length: %d" % int(genome_length), xy=(2,2), xytext=(pi/2, -rmax*0.8), textcoords='data', fontsize=12, horizontalalignment='center', verticalalignment='center')
    ax.annotate("Gene Targets: 6", xy=(2,2), xytext=(pi/2, -rmax*0.95), textcoords='data', fontsize=12, horizontalalignment='center', verticalalignment='center')
    ax.annotate("Status: 6/6 Functional", xy=(2,2), xytext=(3*pi/2, -rmax*0.9), textcoords='data', fontsize=12, horizontalalignment='center', verticalalignment='center')

    # get domains patches
    target_dict = genome.get_targets_from_genome()
    for domain_key in target_dict.keys():
        # create domain 'patch'
        domain = genome.domains[domain_key]
        theta_start = (theta_init + theta_direction * 2.0 * pi * domain.domain_start / genome_length) % (2*pi)
        theta_length = (2.0 * pi * (domain.domain_end - domain.domain_start) / genome_length) % (2*pi)
        ax.bar(theta_start, domain_radius, width=theta_direction*theta_length, color=domain_colours[domain.domain_type], ec='k', alpha=0.75)
        # write domain label
        theta_mid = theta_start + theta_direction * theta_length / 2
        ax.plot(theta_mid, domain_radius, 'o', color='y')
        ax.plot(theta_mid, 0, 'o', color='y')
        plt.text(theta_mid, text_radius, "%s:\n%s" % (domain.domain_type, domain.label), fontsize=10, rotation=(theta_mid*180/pi)+270, horizontalalignment='center',
        verticalalignment='center')

    plt.show()



# TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# initialize genome
from init_genome_camv import init_genome_camv, init_targets_all_domains
complex_concentration = 135000000000
dt = 0.1
pseudo_targets = init_targets_all_domains(complex_concentration)
genome_camv = init_genome_camv(pseudo_targets)
genome_camv.initialize_target_cut_probabilities(dt)
print "plots.."
genome_plot_polar(genome_camv)