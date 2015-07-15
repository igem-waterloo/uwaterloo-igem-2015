"""
Gathering statistics on the Cas9 mutants that were able to bind alternative PAMs

TODO:
- visualization = regions + AA class
- mutations which always co-occur (regionally? or by AA class?)
- mutations which co-occur between NGA and NGC PAMs
- mutations which differ between NGA and NGC PAMs
- sort secondary structure elements by location rather than mutated frequency

"""
import operator
import numpy as np
import matplotlib.pyplot as plt

from cas9_mutants import *


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
    axarr[1].hist(counts_NGC, bins=range(2, 12), histtype='stepfilled', normed=True, align='left', color='#71cce6',
                  label='NGC')
    axarr[1].set_xlabel('Number of Mutations')
    axarr[1].set_ylabel('Probability')
    plt.xlim( 1.5, 10.5 )
    plt.show()


def bar_ss_counts(counts_NGA, counts_NGC):
    """
    Plots a single figure containing two bar charts showing the number of mutations per NGA- or NGA-binding mutant Cas9
    """

    # Sort ss_counts_NGA by number of mutations and ss_counts_NGC to match its order
    counts_NGA = sorted(counts_NGA.items(), key=operator.itemgetter(1))
    sorted_order_NGA = [ss[0] for ss in counts_NGA]
    counts_NGC = sorted(counts_NGC.items(), key=lambda x:sorted_order_NGA.index(x[0]))

    # Bar plot of secondary structure regions containing mutants in each mutant Cas9
    bar_width = 0.45
    plt.bar(np.arange(len(counts_NGA))-bar_width, [ss[1] for ss in counts_NGA], width=bar_width, align='center',
            color='#71cce6', label='NGA')
    print counts_NGC
    plt.bar(range(len(counts_NGC)), [ss[1] for ss in counts_NGC], width=bar_width, align='center',
            color='#333333', label='NGC')
    plt.xticks(range(len(counts_NGA)), [ss[0] for ss in counts_NGA], rotation=45, ha='right')
    plt.title("Mutations in Secondary Structures of Successful Cas9 Mutants")
    plt.xlabel("Secondary Structure")
    plt.ylabel("Mutation Count")
    plt.yscale('log')
    plt.legend()
    plt.show()


def main():

    # Vectors with one entry per mutant Cas9 containing the number of aa mutations in the mutant, separated by PAM
    mutation_counts_NGA = find_mutation_counts('NGA')
    mutation_counts_NGC = find_mutation_counts('NGC')

    # Dictionaries with one key per secondary structure in the PI domain & values that count the number of times the
    # secondary structure was mutated across all mutant Cas9s
    ss_counts_NGA = find_ss_counts('NGA')
    ss_counts_NGC = find_ss_counts('NGC')

    # Plot two histograms & two bar charts
    hist_mutation_counts(mutation_counts_NGA, mutation_counts_NGC)
    bar_ss_counts(ss_counts_NGA, ss_counts_NGC)

if __name__ == '__main__':
    main()
