import scipy
from scipy.stats import kendalltau
from scipy.stats import chisquare


def csv_parser(los):
    """
    Consumes a list of strings containing newline characters and commas, and 
    cleans the strings by removing the newline characters and splitting by 
    commas before returning a list of lists, where each list corresponds to
    the important data from each string, helper for sort_data
    """
    newline_remove = map(lambda (el): el.strip(), los)
    comma_split = map(lambda (el): el.split(','), newline_remove)
    as_ints = [comma_split[0]] + map(lambda (el): map(float, el), comma_split[1:])
    return as_ints


def sort_data(los):
    """
    Consumes a list produced by reading a csv file with titled columns
    corresponding to PyRosetta outputs
    """
    parsed_data = csv_parser(los)
    titles = parsed_data[0]
    just_vals = parsed_data[1:]
    data_transposed = [[just_vals[row][col] for row in xrange(len(just_vals))] for col in xrange(len(titles))]
    sorted_dict = {titles[i]: data_transposed[i] for i in xrange(len(titles))}
    return sorted_dict
        

def calc_tau_and_chi(sorted_data, exp_key):
    """
    Consumes a dictionary of results that is output from sort_data
    and calculates kendall's tau and the chi square statistic value for each
    key with the specified key as the data representing the expected value
    """
    chitau_dict = {}
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    t_exp_array = scipy.array(sorted_data[exp_key])
    c_exp_lst = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for val in sorted_data[exp_key]:
        for i in xrange(1,11):
            if bins[i-1] < val and bins[i] > val:
                c_exp_lst[i-1] += 1
    c_exp_array = scipy.array(c_exp_lst)
    for k in sorted_data.keys():
        k_dict = {}
        t_obs_array = scipy.array(sorted_data[k])
        c_obs_lst = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for val in sorted_data[k]:
            for i in xrange(1,11):
                if bins[i-1] < val and bins[i] > val:
                    c_obs_lst[i-1] += 1
        c_obs_array = scipy.array(c_obs_lst)
        chi_test = chisquare(c_obs_array, f_exp=c_exp_array)
        k_dict['tau'] = kendalltau(t_exp_array, t_obs_array)[0]
        k_dict['chi_sq_val'] = chi_test[0]
        k_dict['chi_sq_p'] = chi_test[1]
        chitau_dict[k] = k_dict
    return chitau_dict


def main(dataname, ref_name, filename):
    """
    Data name is the csv where the data is held, filename is the output file
    name, ref_name is the name of the reference column in the csv
    """
    data = file(dataname, 'r')
    clean_data = sort_data(data.readlines())
    res = calc_tau_and_chi(clean_data, ref_name)
    res_lst = []
    for k in res.keys():
        line = k + ',' + str(res[k]['tau']) + ',' + str(res[k]['chi_sq_val']) \
             + ',' + str(res[k]['chi_sq_p']) + '\n'
        res_lst.append(line)
    res_lst = ['Setup,Tau,ChiSqValue,ChiSqProb\n'] + res_lst
    res_file = file(filename, 'w')
    res_file.writelines(res_lst)
    res_file.close()
