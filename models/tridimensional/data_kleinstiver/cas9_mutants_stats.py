"""
Gathering statistics on the Cas9 mutants that were able to bind alternative PAMs

TODO:
- mutations which always co-occur (regionally? or by AA class?)
- mutations which co-occur between NGA and NGC PAMs
- mutations which differ between NGA and NGC PAMs
"""
import operator
import numpy as np
import matplotlib.pyplot as plt

from cas9_mutants import *
from aa_info import *


def find_mutation_counts(pam):
    """
    Number of mutations that compose each mutant Cas9, separated by PAM
    Returns a vector of the same length as the number of records in mutants_kleinstiver that have 'pam' = pam
    """
    mutation_counts = []
    for mutant in mutants_kleinstiver:
        if mutant['pam'] == pam:
            mutation_counts.append(len(mutant['mutations']))
    return mutation_counts


def find_idx_counts(pam, start=1099, end=1368):
    """
    Cumulative counts of frequency at which amino acid indices are mutated, separated by PAM
    Returns a dict w/ the same keys for indices start-end+1 and values = # of times it was mutated for 'pam' = pam
    """
    idx_counts = {key: 0 for key in range(start, end+1)}
    for mutant in mutants_kleinstiver:
        if mutant['pam'] == pam:
            for mutation in mutant['mutations']:
                if start <= mutation['aa_idx'] <= end:
                    idx_counts[mutation['aa_idx']] += 1
    return idx_counts


def find_ss_counts(pam):
    """
    Cumulative counts of which secondary structures are mutated across mutant Cas9s, separated by PAM
    Returns a dict w/ the same keys as PI_sec_structure and values = # of times it was mutated for 'pam' = pam
    """
    ss_counts = {key: 0 for key in PI_sec_structure}
    for mutant in mutants_kleinstiver:
        if mutant['pam'] == pam:
            for mutation in mutant['mutations']:
                ss_counts[mutation['sec_structure']] += 1
    return ss_counts


def find_aa_changes(pam, change_type, by_idx=True, start=1099, end=1368):
    """
    Interfaces with the aa_info dictionary to classify the amino acid changes across mutant Cas9s, separated by PAM
    if by_idx = True, returns a dict with the keys for each combination of index and aa_info change, otherwise returns
    a dict with keys for every aa_info change observed
    """
    aa_changes = {}
    assert change_type in aa_info.keys(), "change_type is not an in set of aa_info keys: %s" % id
    for mutant in mutants_kleinstiver:
        if mutant['pam'] == pam:
            for mutation in mutant['mutations']:
                aa_info_init = aa_info[change_type][mutation['aa_init']]  # info corresponding to WT nucleotide
                aa_info_mut = aa_info[change_type][mutation['aa_mut']]  # info corresponding to mutated nucleotide

                # key for aa_changes dictionary determined with by_idx setting
                change_key = aa_info_init + 'to' + aa_info_mut
                if by_idx: change_key = str(mutation['aa_idx']) + '_' + change_key

                if change_key in aa_changes:
                    aa_changes[change_key] += 1
                else:
                    aa_changes[change_key] = 1
    return aa_changes


def find_num_pcr_needed():
    """
    Estimate the number of PCR needed to generate each Cas9 mutant from WT Cas9, assuming that WT Cas9 is mutated via
    PCR reactions that can alter 60 nt at a time (or 20 AA).
    Returns: list of ints of the same length as mutants_kleinstiver
    """
    num_pcr_per_mutant = []
    for mutant in mutants_kleinstiver:
        idx_to_pcr_mutate = [m['aa_idx'] for m in mutant['mutations']] # Get indices of all mutations in this mutant
        num_pcr = 0

        # Step inwards from min/max amino acid indices by 60 nt (1 PCR on either side) and remove all indices covered
        while idx_to_pcr_mutate:
            if max(idx_to_pcr_mutate) - min(idx_to_pcr_mutate) > 20: # 2 or more PCR needed
                right_coverage = min(idx_to_pcr_mutate) + 20
                left_coverage = max(idx_to_pcr_mutate) - 20
                idx_to_pcr_mutate = [idx for idx in idx_to_pcr_mutate if right_coverage <= idx <= left_coverage]
                num_pcr += 2
            else:
                num_pcr += 1
                idx_to_pcr_mutate = []  # if the max/min are less than 20 AA apart, only 1 more PCR needed, so end

        num_pcr_per_mutant.append(num_pcr)
    return num_pcr_per_mutant


def hist_mutation_counts(counts_NGA, counts_NGC):
    """
    Plots a single figure containing two histograms showing the number of mutations per NGA- or NGA-binding mutant Cas9
    """
    # Plot histogram of mutations per mutant Cas9
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].hist(counts_NGA, bins=range(2, 12), histtype='stepfilled', normed=True, align='left', color='#71cce6',
                  label='NGA')
    axarr[0].set_title('Mutations per Successful Cas9 Mutant')
    axarr[0].set_ylabel('Probability')
    axarr[0].legend()
    axarr[1].hist(counts_NGC, bins=range(2, 12), histtype='stepfilled', normed=True, align='left', color='#71cce6',
                  label='NGC')
    axarr[1].set_xlabel('Number of Mutations')
    axarr[1].set_ylabel('Probability')
    axarr[1].legend()
    plt.xlim( 1.5, 10.5 )
    plt.show()


def bar_graph_dict(dict_to_plot, plot_title="", xlab="", ylab="", log_scale=False, col="#71cce6", sort_key_list=None,
                   min_count=1):
    """
    Plots a bar graph of the provided dictionary.
    Params:
    dict_to_plot (dict): should have the format {'label': count}
    plot_title (str), xlab (str), ylab (str), log_scale (bool): fairly self-explanatory plot customization
    col (str): colour for the bars
    sort_key_list (list): the keys of the dictionary are assumed to match based its first item and reordered as such
    min_count (int): do not plot items with less than this in the count
    """
    # Sort dictionary & convert to list using custom keys if needed
    if not sort_key_list:
        list_to_plot = sorted(dict_to_plot.items())
    else:
        list_to_plot = sorted(dict_to_plot.items(), key=lambda x: sort_key_list.index(x[0]))

    # Remove list items with less than min_count
    if min_count != 1:
        list_to_plot = [dd for dd in list_to_plot if dd[1] >= min_count]

    # Bar plot of secondary structure regions containing mutants in each mutant Cas9
    bar_width = 0.45
    plt.bar(np.arange(len(list_to_plot)), [dd[1] for dd in list_to_plot], width=bar_width, align='center', color=col)
    plt.xticks(range(len(list_to_plot)), [dd[0] for dd in list_to_plot], rotation=45, ha='right')
    plt.title(plot_title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if log_scale:
        plt.yscale('log')
        plt.ylim(min_count-0.1)  # Show values with just the minimum count
    plt.show()


def main():
    # Vectors with one entry per mutant Cas9 containing the number of aa mutations in the mutant, separated by PAM
    mutation_counts_NGA = find_mutation_counts('NGA')
    mutation_counts_NGC = find_mutation_counts('NGC')

    # Dictionaries with one key per secondary structure in the PI domain & values that count the number of times the
    # secondary structure was mutated across all mutant Cas9s
    ss_counts_NGA = find_ss_counts('NGA')
    ss_counts_NGC = find_ss_counts('NGC')

    # Dictionaries with one key per amino acid PI domain (well, 1097-1364 since those are the values we exported dists
    # for in PyRosetta & values that count the number of times the index was mutated across all mutant Cas9s
    idx_counts_NGA = find_idx_counts('NGA', 1097, 1364)
    idx_counts_NGC = find_idx_counts('NGC', 1097, 1364)

    # Plot two histograms & two bar charts
    hist_mutation_counts(mutation_counts_NGA, mutation_counts_NGC)

    # Plot bar charts of secondary structure, sorted by secondary structure index
    PI_sec_sorted = sorted(PI_sec_structure.items(), key=operator.itemgetter(1))
    PI_sec_sorted = [ss[0] for ss in PI_sec_sorted]
    bar_graph_dict(ss_counts_NGA, log_scale=True, sort_key_list=PI_sec_sorted, ylab="Number of Mutations",
                   plot_title="Mutations Across Secondary Structures for NGA Cas9")
    bar_graph_dict(ss_counts_NGC, log_scale=True, sort_key_list=PI_sec_sorted, ylab="Number of Mutations",
                   plot_title="Mutations Across Secondary Structures for NGC Cas9")

    # Estimate number of PCR reactions needed to convert WT Cas9 into each Cas9 mutant
    num_pcr_needed = find_num_pcr_needed()
    print num_pcr_needed

    # Bar graphs of amino acid category changes, sorted by index
    size_transitions_NGA = find_aa_changes('NGA', 'size', by_idx=False)
    bar_graph_dict(size_transitions_NGA, log_scale=True, ylab="Number of Mutations",
                   plot_title="AA Size Transitions by Index for NGC Cas9", min_count=2)

    size_transitions_NGC = find_aa_changes('NGC', 'size')
    bar_graph_dict(size_transitions_NGA, log_scale=True, ylab="Number of Mutations",
                   plot_title="AA Size Transitions by Index for NGC Cas9", min_count=2)

    # Look at direct amino acid transitions
    aa_transitions_NGA = find_aa_changes('NGA', 'abbreviation')
    bar_graph_dict(aa_transitions_NGA, log_scale=True, ylab="Number of Mutations",
                   plot_title="AA Transitions by Index for NGA Cas9", min_count=2)
    aa_transitions_NGC = find_aa_changes('NGC', 'abbreviation')
    bar_graph_dict(aa_transitions_NGC, log_scale=True, ylab="Number of Mutations",
                   plot_title="AA Transitions by Index for NGC Cas9", min_count=2)


if __name__ == '__main__':
    main()
