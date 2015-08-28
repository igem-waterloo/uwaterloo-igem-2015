import csv
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac

from constants import CSV_HEADER
from utility import int_from_pam_string, pam_string_from_int


def csv_load(fullpath):
    """Loads csv data into memory
    Args:
        fullpath: full path to file (e.g. root/some_dir/some_file.csv)
    Returns:
        tuple: header, and a list of lists created from file.readlines() method
    """
    assert fullpath[-3:] == 'csv'
    with open(fullpath, 'rb') as f:
        reader = csv.reader(f)
        csv_data = []
        for i, row in enumerate(reader):
            if i == 0:
                assert CSV_HEADER[0] == row[0]  # make sure file has header and that it matches expected header start
                csv_header = row
            else:
                csv_data.append(row)
    return csv_header, csv_data


def csv_to_dict(fullpath, keys=['Final DNA']):
    """Loads csv data into memory, then converts the data into dictionary format
    Args:
        fullpath: full path to file (e.g. root/some_dir/some_file.csv)
        keys: list of keys (referencing csv column stats like 'Final DNA' score) to make up the dictionary
    Returns:
        csv in dictionary format where stats reference list of tuples: {statname: (pam, stat)} -- see example
    Example dictionary:
        {'Final DNA':
            ('aaaa', 1234.56),
             ...
            ('tttt', 4321.65)
        }
    """
    csv_header, csv_data = csv_load(fullpath)
    column_index_dict = {key: csv_header.index(key) for key in keys}  # select columns for referencing data
    pam_indices = [i for i, elem in enumerate(csv_header) if 'PAM_' in elem]  # use to concatenate pam columns
    return {key: [(''.join([row[i] for i in pam_indices]),  # concatenate pam
                   float(row[column_index_dict[key]]))      # get stat score corresponding to pam (format as float)
                  for row in csv_data]
            for key in keys}


def sort_csv_data(list_of_tuples, reverse_flag=False):
    """Sort a list of (pam, score) tuples
    Args:
        list_of_tuples: list of tuples of the format [(str, float), ... , (str, float)]
        reverse_flag: [default: False] if True, sort descending instead of ascending
    Returns:
        sorted data in same format
    Notes:
        - sorts by score (second tuple element) in ascending order
    """
    return sorted(list_of_tuples, key=lambda tup: tup[1], reverse=reverse_flag)


def cluster_csv_data(csv_dict, stat_to_cluster='Final DNA'):
    """
    Args:
        csv_dict: dictionary returned from the csv_to_dict() function
        stat_to_cluster: [default: 'Final DNA'] key for csv_dict corresponding to a statistic
    Returns:
        csv data for that statistic in a clustered dictionary format (see example)
    Example:
        key represents the cluster order, in increasing order of sorting
        {1: [(pam, score), ... , (pam, score)]   <-- note: this list is sorted
         2: [(pam, score), ... , (pam, score)],
         ...
         n: [(pam, score), ... , (pam, score)]
        }
    Notes:
        - each list is of arbitrary length
        - n is determined by the clustering algorithm (n <= number of elements in csv data), but can be enforced
    """
    csv_data_as_tuples = csv_dict[stat_to_cluster]
    sorted_data = sort_csv_data(csv_data_as_tuples)
    data_to_cluster = [[int_from_pam_string(pair[0]), pair[1]] for pair in sorted_data]  # convert pams to ints
    cluster_linkage = hac.linkage(data_to_cluster)
    return cluster_linkage


# ====================================================================================
# ====================================================================================
# ====================================================================================
csv_dict = csv_to_dict("Chimera.csv")
print csv_dict
print sort_csv_data(csv_dict['Final DNA'])
linkage = cluster_csv_data(csv_dict)
f = plt.figure()
cluster_dendrogram = hac.dendrogram(linkage)
plt.show()
print "and the linkage...\n"
print linkage
