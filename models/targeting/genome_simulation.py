from init_genome_camv import init_genome_camv


# simulation parameters (time in seconds)
complex_concentration = 135000000000
dt = 0.001
t0 = 0.0
t1 = 3600.0
total_turns = int((t1 - t0) / dt)

# initialize genome to simulate
genome_camv = init_genome_camv(dt=dt, complex_concentration=complex_concentration)


time_sim = t0
for turn in xrange(total_turns):

    # perform simulation actions and collect data
    targets_from_genome = genome_camv.get_targets_from_genome()
    # act on the targets...

    # save turn data (maybe only if stuff happened?)
    # do...

    # update plots if actively showing plots
    # do...

    # increment timer
    time_sim += dt
