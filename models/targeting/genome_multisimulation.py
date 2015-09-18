import datetime
import os

from genome_simulation import genome_simulate


def genome_multisimulate(n):
    multiruns_folder = "runsmulti" + os.sep  # store timestamped runs here
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
    time_folder = current_time + os.sep
    current_run_folder = multiruns_folder + time_folder
    data_folders_basename = "data_"
    dir_list = [multiruns_folder, current_run_folder]
    for dirs in dir_list:
        if not os.path.exists(dirs):
            os.makedirs(dirs)
    for i in xrange(n):
        data_output_path = current_run_folder + data_folders_basename + str(i)
        genome_simulate(flag_plot=False, flag_multirun=True, batch_data_path=data_output_path)
        print "Run %d complete" % i


if __name__ == '__main__':
    genome_multisimulate(1000)
