from numpy import prod
from random import randint

# constants which may be used in other scripts/the simulation itself
# vector of mismatch values with index relative to distance from PAM (Hsu et al. from MIT)
mismatch_decay_values = [0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
                         0.445, 0.508, 0.613, 0.841, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
# average time for SpyCas9 cutting in perfect match case
average_cut_time = 60.0

def prob_concentration(concentration, num_mismatches):
    """Computes the concentration factor of the total cut probability
    Args:
        concentration: intranuclear gRNA-cas9 steady state concentration in nM
        num_mismatches: mismatches between gRNA and target
    Returns:
        float in [0.0, 1.0]
    Notes:
        - the formulas were generated using multiple regression in R with data from Kuscu et al 2014
    """
    if num_mismatches == 0:  # no mismatches case
        return (10**(-0.9935))*(concentration**0.465691)
    else:
        return (10**(0.036775*num_mismatches-1.616842))*(concentration**0.29122)


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
    prob_factor_time = dt/average_cut_time  # assume uniform distribution in time - TODO change to poisson
    mismatch_decay_subset = [mismatch_decay_values[idx] for idx in xrange(len(grna)) if grna[idx] != target[idx]]
    prob_factor_mismatch = prod(mismatch_decay_subset)
    prob_factor_concentration = prob_concentration(concentration, len(mismatch_decay_subset))
    return prob_factor_time * prob_factor_mismatch * prob_factor_concentration

# random insertion function
# takes insertion size and produces random insertion string

DNA_ALPHABET = "acgt"

def nt_rand(insertion_size):
    insertion = ""
    for x in range(insertion_size):
        insertion += DNA_ALPHABET[randint(0,3)]

    return insertion

def indel():
    del_left = randint(0,4)
    del_right = randint(0,4)
    insert = randint(0,8)

    return [del_left, del_right, insert]
