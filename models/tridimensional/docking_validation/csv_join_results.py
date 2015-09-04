import argparse
import csv
import datetime
import fnmatch
import os.path

from collections import defaultdict
from constants import CSV_HEADER, PAM_TOOLS


def locate(pattern, root=os.curdir):
    """Locate all files matching supplied filename pattern in and below supplied root directory.
    From http://code.activestate.com/recipes/499305-locating-files-throughout-a-directory-tree/
    """
    print pattern
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


if __name__ == '__main__':
    # create parser and parse arguments
    parser = argparse.ArgumentParser(description='Join CSV results from a directory into a single file')
    parser.add_argument('--csv_dir', metavar='D', default='', type=str,
                        help='path to directory containing CSVs to merge (default: ".")')
    parser.add_argument('-o', '--output_dir', metavar='D', type=str,
                        help='path to output directory for merged CSV')
    parser.add_argument('--merge_stat', metavar='S', default="Final DNA", type=str,
                        help='string of CSV column to merge, matching header (default: "Final DNA")')
    parser.add_argument('--alt_tool', metavar='S', nargs='?', const='3DNA', default='Chimera',
                        type=str, help='[switch] use 3DNA csvs instead of Chimera (default: Chimera)')
    args = parser.parse_args()

    assert args.merge_stat in CSV_HEADER, "--merge_stat %r is not a recognized CSV header column" % args.merge_stat
    assert os.path.isdir(args.csv_dir), "--csv_dir %r is not a directory, check provided option" % args.csv_dir
    assert args.alt_tool in PAM_TOOLS, "--alt_tool not in list of accepted PAM tools"

    # Setup output path for scoring
    if args.output_dir is not None:
        csv_output_dir = args.output_dir
    else:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %I.%M.%S%p")
        csv_output_dir = os.path.join("results", timestamp)
    os.makedirs(csv_output_dir)

    # For each CSV, the column containing the merge_stat will be associated in a dictionary with the PAM
    merge_stat_idx = CSV_HEADER.index(args.merge_stat)
    joined_colnames = ['pam']
    csvs_to_join = []

    # Find all CSVs in directory
    print "Locating CSVs in "+args.csv_dir+"..."
    for csv_to_join in locate("*"+args.alt_tool+"*.csv", args.csv_dir):

        # Assume the directory containing the csv has a descriptive name, use as column name
        joined_colnames.append(os.path.basename(os.path.dirname(csv_to_join)))
        dict_to_append = {}

        # Read every CSV for the column of interest
        with open(csv_to_join) as c:
            r = csv.reader(c, delimiter=',')
            for row in r:
                pam_key = ''.join(row[0:4])
                merge_stat_val = row[merge_stat_idx]
                dict_to_append[pam_key] = merge_stat_val

        csvs_to_join.append(dict(dict_to_append))

    # Join the list of dictionaries on chosen column
    joined_csv = defaultdict(list)
    for d in csvs_to_join:
        for key, value in d.iteritems():
            joined_csv[key].append(value)

    # Print to a new CSV
    with open(os.path.join(csv_output_dir, 'joined_csv_results.csv'), 'wb') as f:
        w = csv.writer(f)
        w.writerow(joined_colnames)  # removes trailing comma
        for key, value in joined_csv.iteritems():
            w.writerow([key] + value)
