import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.cluster.hierarchy as hac
import scipy.spatial.distance as scidist
import shutil

from constants import CSV_HEADER
from utility import int_from_pam_string, pam_string_from_int


def csv_load(fullpath):
    """Loads csv data into memory
    Args:
        fullpath: full path to file (e.g. root/some_dir/some_file.csv)
    Returns:
        tuple: header, and a list of lists created from file.readlines() method
    """
    assert fullpath[-4:] == '.csv'
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
        csv in dictionary format where stats reference dictionary: {statname: {pam: stat value} } -- see example
    Example dictionary:
        {'Final DNA':
            {'aaaa': 1234.56,
             ...
            'tttt': 4321.65}
        }
    """
    csv_header, csv_data = csv_load(fullpath)
    column_index_dict = {key: csv_header.index(key) for key in keys}  # select columns for referencing data
    pam_indices = [i for i, elem in enumerate(csv_header) if 'PAM_' in elem]  # use to concatenate pam columns
    print pam_indices
    csv_dict = {key: {''.join([row[i] for i in pam_indices]): float(row[column_index_dict[key]])
                      for row in csv_data}  # concatenate pam / get stat score corresponding to pam (format as float)
                for key in keys}  # do this for all the desired statistics
    csv_dict['header'] = csv_header
    return csv_dict


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
    csv_data_as_keyvalue = csv_dict[stat_to_cluster]
    data_to_cluster = [[csv_data_as_keyvalue[key]] for key in csv_data_as_keyvalue.keys()]
    cluster_linkage = hac.linkage(data_to_cluster, method='single', metric='euclidean')
    return cluster_linkage


def plot_cluster_dendrogram(cluster_linkage, length_pam, threshold='default'):
    """Dendrograms are representations of heirarchical clusters
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    leaf_label_map = lambda x: pam_string_from_int(x, length_pam)
    plt.figure()
    dendrodata = hac.dendrogram(cluster_linkage, color_threshold=threshold, leaf_label_func=leaf_label_map, leaf_rotation=45.0, leaf_font_size=8)
    # segment to plot distances between clusters
    """
    for i, d in zip(dendrodata['icoord'], dendrodata['dcoord']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'ro')
        plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                     textcoords='offset points',
                     va='top', ha='center')
    """
    plt.show()


def cluster_csv_data(csv_dict, stat_to_cluster='Final DNA', plot_dendrogram_flag=True):
    """Clusters linkage object by applying a threshold to get a flat clustering
    Args:
        csv_dict: dictionary returned from the csv_to_dict() function
        stat_to_cluster: [default: 'Final DNA'] key for csv_dict corresponding to a statistic
        plot_dendrogram_flag: plot dendrogram if True
    Returns:
        csv data for that statistic in a clustered dictionary format (see example)
    Example of returned dictionary:
        {pam:
            {'stat_value': float,  <-- data value that's been clustered
             'stat_cluster': int,  <-- cluster rank (1 to n)
             'stat_cluster_centroid': float},  <-- average value of associated cluster rank
         ... }  <-- for all pams
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    # prepare data and compute distances
    csv_data_as_keyvalue = csv_dict[stat_to_cluster]
    length_pam = len(csv_data_as_keyvalue.keys()[0])
    print csv_data_as_keyvalue.keys()
    print length_pam
    length_data = len(csv_data_as_keyvalue.keys())
    data_to_cluster = [[csv_data_as_keyvalue[key]] for key in csv_data_as_keyvalue.keys()]  # ignore pams, keep order
    pair_dists = scidist.pdist(data_to_cluster, metric='euclidean')
    # determine cluster membership
    linkage = get_cluster_linkage(csv_dict, stat_to_cluster=stat_to_cluster)
    threshold = 0.5 * np.std(pair_dists)
    cluster_membership_array = hac.fcluster(linkage, threshold, criterion='distance')
    print threshold
    print cluster_membership_array
    print set(cluster_membership_array)
    # revise cluster membership so that 'best' means cluster 1 instead of the last cluster (lower energy is better)
    cluster_low, cluster_high = np.min(cluster_membership_array), np.max(cluster_membership_array)
    cluster_sum = cluster_low + cluster_high
    transform = lambda val: cluster_sum - val
    for i, elem in enumerate(cluster_membership_array):
        cluster_membership_array[i] = transform(elem)
    # assign cluster membership
    clustered_data = {}
    for i, key in enumerate(csv_data_as_keyvalue.keys()):
        clustered_data[key] = {'stat_value': csv_data_as_keyvalue[key],
                               'stat_cluster': cluster_membership_array[i],  # TODO CHECK ORDER?
                               'stat_cluster_centroid': None,  # TODO CHECK ORDER? AND IMPLEMENT
                               'pam_sanity_check': pam_string_from_int(i, length_pam)}  # TODO FIX AND REMOVE
    # conditionally plot dendrogram
    if plot_dendrogram_flag:
        plot_cluster_dendrogram(linkage, length_pam, threshold=threshold)
    return clustered_data


def write_clustered_csv(fullpath_input, dir_output=None, stats_to_cluster=['Final DNA']):
    """Clusters specific data from an input csv and writes a new csv with appended clustering information
    Args:
        fullpath_input: full path to the input csv
        dir_output: directory where the output csv will be placed
        stats_to_cluster: list of stats to cluster
    Returns:
        full path to new csv with appended clustering information
    """
    # IO preparation
    assert fullpath_input[-4:] == '.csv'
    dirpath, filename_input = os.path.split(fullpath_input)
    filename_output = filename_input[:-4] + '_clustered.csv'
    if dir_output is None:
        fullpath_output = os.path.join(dirpath, filename_output)
    else:
        fullpath_output = os.path.join(dir_output, filename_output)

    # load data for clustering
    csv_dict = csv_to_dict(fullpath_input, keys=stats_to_cluster)
    csv_header = csv_dict['header']
    pam_indices = [i for i, elem in enumerate(csv_header) if 'PAM_' in elem]  # use to concatenate pam columns

    # cluster each stat separately and track header changes
    cluster_dict = {}
    csv_cluster_header = []
    for stat in stats_to_cluster:
        cluster_dict[stat] = cluster_csv_data(csv_dict, stat_to_cluster=stat, plot_dendrogram_flag=False)
        csv_cluster_header.append('%s cluster' % stat)
        # csv_cluster_header.append('%s cluster centroid' % stat)  # TODO implement cluster centroid

    # write clustered data to csv
    data_to_append = ['stat_cluster']  # , 'stat_cluster_centroid']  # TODO implement cluster centroid
    with open(fullpath_input, 'r') as csvin:
        reader = csv.reader(csvin)
        with open(fullpath_output, 'wb') as csvout:
            writer = csv.writer(csvout)
            for i, row in enumerate(reader):
                if i == 0:
                    writer.writerow(row + csv_cluster_header)
                else:
                    pam = ''.join([row[i] for i in pam_indices])
                    cluster_data_to_append = []
                    for stat in stats_to_cluster:
                        cluster_data_to_append += [cluster_dict[stat][pam][key] for key in data_to_append]
                    writer.writerow(row + cluster_data_to_append)

    return fullpath_output


if __name__ == '__main__':
    # TODO sys.argv
    #csv_dict = csv_to_dict("Chimera.csv")
    csv_dict = csv_to_dict("Chimera_64.csv")
    clustered_data = cluster_csv_data(csv_dict, plot_dendrogram_flag=True)
    print clustered_data

    #write_clustered_csv("Chimera.csv")
    write_clustered_csv("Chimera_64.csv")
