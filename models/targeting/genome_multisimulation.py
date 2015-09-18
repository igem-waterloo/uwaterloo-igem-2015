import datetime
import os

from genome_csv import multirun_gene_state_compile_to_dict, multirun_gene_state_compile_to_csv
from genome_plot import plot_state_multirun
from genome_simulation import genome_simulate


def genome_multisimulate(n):

    # io setup
    multiruns_folder = "runsmulti" + os.sep  # store timestamped runs here
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
    time_folder = current_time + os.sep
    current_run_folder = multiruns_folder + time_folder
    data_folders_basename = "data_"
    dir_list = [multiruns_folder, current_run_folder]
    for dirs in dir_list:
        if not os.path.exists(dirs):
            os.makedirs(dirs)

    # perform n simulations
    for i in xrange(n):
        data_output_path = current_run_folder + data_folders_basename + str(i)
        genome_simulate(flag_plot=False, flag_multirun=True, batch_data_path=data_output_path)
        if i % (n/10) == 0:
            print "Run %d/%d complete" % (i, n)  # plot progress every 10% or so

    # compile n simulations
    state_totals_gene = multirun_gene_state_compile_to_dict(current_run_folder)
    multirun_gene_state_compile_to_csv(current_run_folder, "state_totals_gene.csv")
    plot_state_multirun(state_totals_gene, "gene", labels_to_plot=["gene_P6"],
                        output_path=os.path.join(current_run_folder, "states_gene_totals.png"), flag_show=False)

    return


if __name__ == '__main__':
    n = 100
    genome_multisimulate(n)
