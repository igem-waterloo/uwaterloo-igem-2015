# probability of cleavage function
# Some comments:
# 1. the cut probability vector (mmval), is from a paper by Hsu et al, as
#    recorded on the MIT CRISPR edu page
# 2. the affect of dt is accounted for by assuming uniform cutting probability
#    over time with a set time to bond and cut of 60 s, this is very arbitrary
#    and based on only spurious mentions in a paper of time, if anyone has
#    better implementation ideas add them
# 3. in combining the efficiencies, I assumed independence of cut probability
#    affects, this is once again arbitrary, and in fact represents a 
#    boundary case to the true interaction, likely the true affects are more
#    highly negatively dependent
# 4. the relationship to concentration was determined using multiple regression
#    on data from at 2014 paper by Kuscu, here they determined the proportion
#    bonded given different concentrations and numbers of matching base pairs,
#    this is a best case, however, as they only did variations on the matches
#    furthest from the PAM site, which are very non-specific
#    this function is separated into the case without mismatches and the case
#    with n mismatches because the case without mismatches is individually
#    important and better fit by a particular level curve

import math

def cut_prob(gRNA, target, dt, concentration):
    mmval = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.841,
             0.732,0.828,0.615,0.804,0.685,0.583] # vector of mismatch values
                                                  # from Hsu et al. from MIT
                                                  # CRISPR site
    
    bond_cut_t = 60 # couldn't find anything other than 0-3 min for time
                    # so plugged in 60 s as a possible value, arbitrary
    dt_prob = dt/bond_cut_t # simple relative length of dt to determine this
                            # probability factor
    
    mmind = []
    for i in range(23): # checking all positions for mismatches
        if gRNA[i] != target[i]:
            mmind.append(i)
    
    rel_prob = []
    for ind in mmind: # extract relevant relative efficiencies from mmval
        rel_prob.append(mmval[ind])
        
    tot_prob = 1
    for p in rel_prob: # collecting probabilities assuming independence(weak)
        tot_prob = tot_prob*p
    
    def conc_prob(conc,n): # here in nM, data from Kuscu et al 2014
        if n == 0: # case if there are no mismatches
            return (10**(-0.9935))*(conc**0.465691)
        else:   
            return (10**(0.036775*n-1.616842))*(conc**0.29122)
        # the formula above was generated using multiple regression in R
    
    conc_p = conc_prob(concentration, len(mmind)) # get concentration prob for
                                                  # particular situation
    
    return tot_prob*conc_p*dt_prob # return product of probs