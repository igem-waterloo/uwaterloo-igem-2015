import math
from numpy import prod
from numpy.random import randint, choice


# constants which may be used in other scripts/the simulation itself
DNA_ALPHABET = "acgt"
# vector of mismatch values with index relative to distance from PAM (Hsu et al. from MIT)
mismatch_decay_values = [0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
                         0.445, 0.508, 0.613, 0.841, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
# average time for SpyCas9 cutting in perfect match case
average_cut_time = 60.0
# distribution of indel size from CRISPResso
indel_sizes = range(-20, 9)
indel_probs = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
               0.001, 0.001, 0.003, 0.003, 0.001, 0.023, 0.001, 0.005, 0.007, 
               0.016, 0.01, 0.06, 0.84, 0.01, 0.005, 0.001, 0.001, 0.001,
               0.001, 0.001, 0.001]
# distribution of insertions from CRISPResso
insertion_sizes = range(0,9)
insertion_probs = [0.43, 0.053, 0.01, 0.035, 0.42, 0.045, 0.005, 0.001, 0.001]

def prob_concentration(concentration, num_mismatches):
    """Computes the concentration factor of the total cut probability
    Args:
        concentration: intranuclear gRNA-cas9 steady state concentration in nM
        num_mismatches: mismatches between gRNA and target
    Returns:
        float in [0.0, 1.0]
    Notes:
        - the formulas were generated using regression and fitting 
        in Excel with data from Wu et al 2014
    """
    glob_scale = 2/math.pi # ensure it's between 0 and 1
    mismatch_scale = 354.3859 - 15.8762*num_mismatches # generated using regression in Excel
    return (glob_scale*math.atan(concentration/mismatch_scale))
    #Holding on to old R regression formulas for now in case the new suggested
    #formula is not acceptable, and for easy reference
    #if num_mismatches == 0:  # no mismatches case
    #    return (10**(-0.9935))*(concentration**0.465691)
    #else:
    #    return (10**(0.036775*num_mismatches-1.616842))*(concentration**0.29122)


def prob_cut(grna, target, concentration, dt):
    """Probability of cleavage function
    Args:
        grna: 20 long string of nucleotides
        target: 20 long string of nucleotides
        concentration: intranuclear gRNA-cas9 steady state concentration in nM
        dt: simulation timestep
    Returns:
        float in [0.0, 1.0] representing the probability of cutting within the window of time dt
    Notes:
        - use Hsu et al. (from MIT) data for the decay probability of cutting with mismatches
        - the affect of dt is accounted for by assuming uniform cutting probability over time with a set time to
          bond and cut of 60 s, this is very arbitrary and based on only spurious mentions in a paper of time,
          if anyone has better implementation ideas add them
        - in combining the efficiencies, I assumed independence of cut probability effects, this is once again
          arbitrary, and in fact represents a boundary case to the true interaction, likely the true affects
          are more highly negatively dependent
        - the relationship to concentration was determined using multiple regression on data from at
          2014 paper by Kuscu, here they determined the proportion bonded given different concentrations
          and numbers of matching base pairs, this is a best case, however, as they only did variations on the
          matches furthest from the PAM site, which are very non-specific this function is separated into
          the case without mismatches and the case with n mismatches because the case without mismatches
          is individually important and better fit by a particular level curve
    """
    assert len(grna) == len(target)
    prob_factor_time = 1 - math.exp(-dt/average_cut_time) 
    mismatch_decay_subset = [mismatch_decay_values[idx] for idx in xrange(len(grna)) if grna[idx] != target[idx]]
    prob_factor_mismatch = prod(mismatch_decay_subset)
    prob_factor_concentration = prob_concentration(concentration, len(mismatch_decay_subset))
    return prob_factor_time * prob_factor_mismatch * prob_factor_concentration


def nt_rand(insertion_size):
    """Takes insertion size and produces random insertion string
    """
    insertion = ""
    for x in xrange(insertion_size):
        insertion += DNA_ALPHABET[randint(0, 3)]
    return insertion

def indel():
    """
    Compute the indel characteristics for use in other functions, outputs
    a list where the elements represent the size of deletions to the left and 
    right and the size of insertion in that order
    To ensure the distribution is true to the CRISPResso data the size of indel
    is first chosen according to the CRISPResso distribution, after which the
    insertion size is chosen according to the CRISPResso distribution. Following
    this the deletions to the left and right are chosen such that the selected
    indel size equals the insertion minus the left and right deletions and such
    that deletions to the left are likely to be larger than those to the right,
    once again to be consistent with CRISPResso data
    """
    size = choice(indel_sizes, 1, p = indel_probs)[0]
    insert = choice(insertion_sizes, 1, p = insertion_probs)[0]
    total_del = abs(size - insert)
    if total_del != 0:
        del_left = randint(math.floor(total_del/4), total_del)
        del_right = total_del - del_left
    else:
        del_left = 0
        del_right = 0
    return [del_left, del_right, insert]
