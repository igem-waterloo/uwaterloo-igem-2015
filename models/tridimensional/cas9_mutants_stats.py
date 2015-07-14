"""
Gathering statistics on the Cas9 mutants that were able to bind alternative PAMs

TODO:
- visualization = regions + AA class
- mutations which always co-occur (regionally? or by AA class?)
- mutations which co-occur between NGA and NGC PAMs
- mutations which differ between NGA and NGC PAMs

"""
import operator
import matplotlib.pyplot as plt

from cas9_mutants import *

# Count number of mutations that compose each mutant Cas9, separated by PAM
mutation_counts_NGA = []
mutation_counts_NGC = []

# Cumulative counts of which secondary structures are mutated across mutant Cas9s, separated by PAM
ss_counts_NGA = {key: 0 for key in PI_sec_structure}
ss_counts_NGC = {key: 0 for key in PI_sec_structure}

for mutant in mutants_kleinstiver:
    if mutant['pam'] == 'NGA':
        mutation_counts_NGA.append(len(mutant['mutations']))
        for mutation in mutant['mutations']:
            ss_counts_NGA[mutation['sec_structure']] += 1
    elif mutant['pam'] == 'NGC':
        mutation_counts_NGC.append(len(mutant['mutations']))
        for mutation in mutant['mutations']:
            ss_counts_NGC[mutation['sec_structure']] += 1

# Plot histogram of mutations per mutant Cas9
plt.hist(mutation_counts_NGA, bins = range(2, 12), histtype='stepfilled', normed=True,
         align='left', color='b', label='NGA')
plt.hist(mutation_counts_NGC, bins = range(2, 12), histtype='stepfilled', normed=True,
         align='left', color='r', alpha=0.5, label='NGC')
plt.title("Mutations per Successful Cas9 Mutant")
plt.xlabel("Number of Mutations")
plt.ylabel("Probability")
plt.xlim( 1.5, 10.5 )
plt.legend()
plt.show()

# Sort ss_counts_NGA by number of mutations and ss_counts_NGC to match its order
ss_counts_NGA = sorted(ss_counts_NGA.items(), key=operator.itemgetter(1))
ss_sorted_order_NGA = [ss[0] for ss in ss_counts_NGA]
ss_counts_NGC = sorted(ss_counts_NGC.items(), key=lambda x:ss_sorted_order_NGA.index(x[0]))

# Bar plot of secondary structure regions containing mutants in each mutant Cas9
plt.bar(range(len(ss_counts_NGA)), [ss[1] for ss in ss_counts_NGA], align='center', color='b', label='NGA')
plt.bar(range(len(ss_counts_NGC)), [ss[1] for ss in ss_counts_NGC], align='center', color='r', alpha=0.5, label='NGA')
plt.xticks(range(len(ss_counts_NGA)), [ss[0] for ss in ss_counts_NGA], rotation=45, ha='right')
plt.title("Mutations in Secondary Structures of Successful Cas9 Mutants")
plt.xlabel("Secondary Structure")
plt.ylabel("Mutation Count")
plt.legend()
plt.show()
