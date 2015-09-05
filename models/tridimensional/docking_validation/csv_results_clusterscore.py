import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.cluster.hierarchy as hac
import scipy.spatial.distance as scidist

from csv_results import csv_header_bugfix
from constants import CSV_HEADER, DNA_ALPHABET
from utility import safe_mkdir


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
    # legacy header bug fix check (64pam run but header is for 256pam)
    if csv_header[3] == 'PAM_4' and csv_data[0][3] not in DNA_ALPHABET:
        csv_header.pop(3)
        csv_header_bugfix(fullpath)
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
            'tttt': 4321.65}}
    """
    csv_header, csv_data = csv_load(fullpath)
    column_index_dict = {key: csv_header.index(key) for key in keys}  # select columns for referencing data
    pam_indices = [i for i, elem in enumerate(csv_header) if 'PAM_' in elem]  # use to concatenate pam columns
    csv_dict = {}
    for key in keys:
        csv_dict[key] = {}
        for row in csv_data:
            pam = ''.join([row[i] for i in pam_indices])  # concatenate pam
            csv_dict[key][pam] = float(row[column_index_dict[key]])  # get pam's stat value
    csv_dict['header'] = csv_header
    return csv_dict


def get_cluster_linkage(data_to_cluster):
    """Gets a linkage object representing heirarchical cluster options defined by distance thresholds
    Args:
        data_to_cluster: vectorized data from dictionary returned from the csv_to_dict() function
        stat_to_cluster: [default: 'Final DNA'] key for csv_dict corresponding to a statistic
    Returns:
        linkage object
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    cluster_linkage = hac.linkage(data_to_cluster, method='single', metric='euclidean')
    return cluster_linkage


def plot_cluster_dendrogram(cluster_linkage, keylist, threshold='default'):
    """Dendrograms are representations of heirarchical clusters
    See documentation:
        http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    """
    leaf_label_map = lambda x: keylist[x]
    plt.figure()
    dendrodata = hac.dendrogram(cluster_linkage, color_threshold=threshold, leaf_label_func=leaf_label_map,
                                leaf_rotation=45.0, leaf_font_size=8)
    # TODO find a way to make the extra information below less cramped when plotting
    """# segment to plot distances between clusters
    for i, d in zip(dendrodata['icoord'], dendrodata['dcoord']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'ro')
        plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                     textcoords='offset points',
                     va='top', ha='center')"""
    if threshold != 'default':
        plt.axhline(threshold, color='k', linestyle='--', label='threshold')
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
    data_to_cluster = [[csv_data_as_keyvalue[key]] for key in csv_data_as_keyvalue.keys()]  # ignore pams, keep order
    pair_dists = scidist.pdist(data_to_cluster, metric='euclidean')

    # determine cluster membership
    linkage = get_cluster_linkage(data_to_cluster)
    threshold = 0.5 * np.std(pair_dists)
    cluster_membership_array = hac.fcluster(linkage, threshold, criterion='distance')

    # compute cluster centroid dictionary
    lens_dict = {cluster_idx: 0 for cluster_idx in set(cluster_membership_array)}
    centroid_dict = {cluster_idx: 0.0 for cluster_idx in set(cluster_membership_array)}
    for i, cluster_idx in enumerate(cluster_membership_array):
        centroid_dict[cluster_idx] += csv_data_as_keyvalue.values()[i]
        lens_dict[cluster_idx] += 1
    for cluster_idx in centroid_dict.keys():
        centroid_dict[cluster_idx] = centroid_dict[cluster_idx] / lens_dict[cluster_idx]  # take average

    # revise cluster membership so that 'best' means cluster 1 instead of the last cluster (lower energy is better)
    order_map = {cluster_idx: 0 for cluster_idx in set(cluster_membership_array)}
    centroid_dict_reversed = {centroid_dict[cluster_idx]: cluster_idx for cluster_idx in centroid_dict.keys()}
    for cluster_idx in set(cluster_membership_array):
        min_centroid = np.min(centroid_dict_reversed.keys())
        order_map[cluster_idx] = centroid_dict_reversed[min_centroid]
        del centroid_dict_reversed[min_centroid]
    order_map = {order_map[key]: key for key in order_map.keys()}  # invert order map so it functions as intended
    transform = lambda cluster_idx: order_map[cluster_idx]
    for i, elem in enumerate(cluster_membership_array):  # transform cluster_membership_array references
        cluster_membership_array[i] = transform(elem)
    centroid_dict = {order_map[cluster_idx]: centroid_dict[cluster_idx] for cluster_idx in centroid_dict.keys()}

    # assign cluster membership
    clustered_data = {}
    for i, key in enumerate(csv_data_as_keyvalue.keys()):
        clustered_data[key] = {'stat_value': csv_data_as_keyvalue[key],
                               'stat_cluster': cluster_membership_array[i],
                               'stat_cluster_centroid': centroid_dict[cluster_membership_array[i]]}

    # conditionally plot dendrogram
    if plot_dendrogram_flag:
        plot_cluster_dendrogram(linkage, csv_data_as_keyvalue.keys(), threshold=threshold)

    return clustered_data


def write_clustered_csv(fullpath_input, dir_output=None, stats_to_cluster=['Final DNA'], plot_dendrogram_flag=False):
    """Clusters specific data from an input csv and writes a new csv with appended clustering information
    Args:
        fullpath_input: full path to the input csv
        dir_output: directory where the output csv will be placed
        stats_to_cluster: list of stats to cluster
        plot_dendrogram_flag: selectively plot the dendrogram of clusters
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
        safe_mkdir(dir_output)
        fullpath_output = os.path.join(dir_output, filename_output)

    # load data for clustering
    csv_dict = csv_to_dict(fullpath_input, keys=stats_to_cluster)
    csv_header = csv_dict['header']
    pam_indices = [i for i, elem in enumerate(csv_header) if 'PAM_' in elem]  # use to concatenate pam columns

    # cluster each stat separately and track header changes
    cluster_dict = {}
    csv_cluster_header = []
    for stat in stats_to_cluster:
        cluster_dict[stat] = cluster_csv_data(csv_dict, stat_to_cluster=stat, plot_dendrogram_flag=plot_dendrogram_flag)
        csv_cluster_header.append('%s cluster' % stat)
        csv_cluster_header.append('%s cluster centroid' % stat)

    # write clustered data to csv
    data_to_append = ['stat_cluster', 'stat_cluster_centroid']
    with open(fullpath_input, 'r') as csvin:
        reader = csv.reader(csvin)
        with open(fullpath_output, 'wb') as csvout:
            writer = csv.writer(csvout)
            for i, row in enumerate(reader):
                if i == 0:
                    writer.writerow(csv_header + csv_cluster_header)
                else:
                    pam = ''.join([row[i] for i in pam_indices])
                    cluster_data_to_append = []
                    for stat in stats_to_cluster:
                        cluster_data_to_append += [cluster_dict[stat][pam][key] for key in data_to_append]
                    writer.writerow(row + cluster_data_to_append)

    print "Finished writing cluster results to %s" % fullpath_output
    return fullpath_output


if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(description='Cluster given csv data into a new clustered csv.')
    parser.add_argument('--path_input', metavar='C', type=str, help='directory of input csv file')
    parser.add_argument('--dir_output', metavar='S', nargs='?', default=None,
                        type=str, help='directory to place output csv (default: same as input)')
    parser.add_argument('--plot_on', metavar='F', nargs='?', const=True, default=False,
                        type=str, help='[switch] plot cluster dendrogram (default: no plot)')
    args = parser.parse_args()
    # write to csv
    write_clustered_csv(args.path_input, dir_output=args.dir_output, plot_dendrogram_flag=args.plot_on)
