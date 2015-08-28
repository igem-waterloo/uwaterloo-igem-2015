import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as hac
import scipy.spatial.distance as scidist

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


def sort_tuples_by_idx(list_of_tuples, tuple_idx=1, reverse_flag=False):
    """Sort a list of (pam, score) tuples
    Args:
        list_of_tuples: list of tuples of the format [(str, float), ... , (str, float)]
        tuple_idx: [default: 1] tuple index which defines sorting
        reverse_flag: [default: False] if True, sort descending instead of ascending
    Returns:
        sorted data in same format
    Notes:
        - sorts by score (second tuple element) in ascending order
    """
    return sorted(list_of_tuples, key=lambda tup: tup[tuple_idx], reverse=reverse_flag)


def get_cluster_linkage(csv_dict, stat_to_cluster='Final DNA'):
    """Gets a linkage object representing heirarchical cluster options defined by distance thresholds
    Args:
        csv_dict: dictionary returned from the csv_to_dict() function
        stat_to_cluster: [default: 'Final DNA'] key for csv_dict corresponding to a statistic
    Returns:
        linkage object
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    csv_data_as_tuples = csv_dict[stat_to_cluster]
    data_to_cluster = [[int_from_pam_string(pair[0]), pair[1]] for pair in csv_data_as_tuples]  # convert pams to ints
    #cluster_linkage = hac.linkage(data_to_cluster, method='single', metric='euclidean')
    cluster_linkage = hac.linkage(data_to_cluster, method='complete', metric='euclidean')
    return cluster_linkage


def plot_cluster_dendrogram(cluster_linkage, length_pam, threshold='default'):
    """Dendrograms are representations of heirarchical clusters
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    leaf_label_map = lambda x: pam_string_from_int(x, length_pam)
    plt.figure()
    hac.dendrogram(cluster_linkage, color_threshold=threshold, leaf_label_func=leaf_label_map, leaf_rotation=45.0, leaf_font_size=8)
    plt.show()


def cluster_csv_data(csv_dict, stat_to_cluster='Final DNA', plot_dendrogram_flag=True):
    """Clusters linkage object by applying a threshold to get a flat clustering
    Args:
        csv_dict: dictionary returned from the csv_to_dict() function
        stat_to_cluster: [default: 'Final DNA'] key for csv_dict corresponding to a statistic
        plot_dendrogram_flag: plot dendrogram if True
    Returns:
        csv data for that statistic in a clustered dictionary format (see example)
    Example:
        key represents the cluster order, in increasing order of sorting
        {1: [(pam, score), ... , (pam, score)]   <-- note: this list is sorted
         2: [(pam, score), ... , (pam, score)],
         ...
         n: [(pam, score), ... , (pam, score)]
        }
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    # prepare data and compute distances
    csv_data_as_tuples = csv_dict[stat_to_cluster]
    length_pam = len(csv_data_as_tuples[0][0])
    length_data = len(csv_data_as_tuples)
    data_to_cluster = [[int_from_pam_string(pair[0]), pair[1]] for pair in csv_data_as_tuples]  # convert pams to ints
    pair_dists = scidist.pdist(data_to_cluster, metric='euclidean')
    print pair_dists
    # determine cluster membership
    linkage = get_cluster_linkage(csv_dict, stat_to_cluster=stat_to_cluster)
    threshold = 0.5 * np.std(pair_dists)
    cluster_membership_array = hac.fcluster(linkage, threshold, criterion='distance')
    print threshold
    print cluster_membership_array
    print set(cluster_membership_array)
    # assign cluster membership
    clustered_data = [0] * length_data
    for i, pair in enumerate(data_to_cluster):
        clustered_data[i] = (pam_string_from_int(pair[0], length_pam), pair[1], cluster_membership_array[i])
    # conditionally plot dendrogram
    if plot_dendrogram_flag:
        plot_cluster_dendrogram(linkage, length_pam, threshold=threshold)
    return clustered_data


# ====================================================================================
# ====================================================================================
# ====================================================================================
csv_dict = csv_to_dict("Chimera.csv")
clustered_data = cluster_csv_data(csv_dict, plot_dendrogram_flag=True)
print clustered_data
print sort_tuples_by_idx(clustered_data, tuple_idx=2)
