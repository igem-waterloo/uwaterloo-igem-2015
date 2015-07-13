import argparse
import csv
import os


# constants
header = ['PAM_1', 'PAM_2', 'PAM_3', 'Tool', 'Init FA', 'Init DNA', 'Final FA', 'Final DNA']
categories = ['3DNA', 'Chimera']  # tools used to create 64 pam variants


def get_score_info(score_dir, score_filename):
    """Extracts pam/tool info and scores as a list from a score file
    Args:
        score_dir: path to score file directory
        score_file: specific score file name (including .txt)
    Returns:
        list of scores (as strings)
    Example:
        Example score file: 'agg_3DNA.txt'
            Final FA Score:  2465.219
            Init DNA Score:   305.504
            Final FA score:  2044.202
            Final DNA score: -207.921
            Time taken:       857.995
        Expected return: ['a', 'g', 'g', '3DNA', '2465.219', ..., '-207.921']
    """
    # get textfile classifiers (PAM and category)
    pam, pam_tool_txt = score_filename.split('_')
    pam_tool = pam_tool_txt.split('.')[0]
    assert pam_tool in categories
    # read contents of textfile
    with open(score_dir + score_filename) as f:
        file_lines = f.readlines()
        assert len(file_lines) == 6
    # format row
    scores = [file_lines[i].split(" ")[-1][:-1] for i in xrange(4)]
    return [pam[0], pam[1], pam[2], pam_tool] + scores


def results_to_csv(score_file_directory):
    """Appends data from each result file to an appropriate csv
    Args:
        score_file_directory: path to score files
    Returns:
        None
    Notes:
        - currently creates / looks for a csv for each PAM tool
    """
    # csv prep
    csv_dict = {elem: csv.writer(open('%s.csv' % elem, 'a'), lineterminator='\n') for elem in categories}
    for elem in categories:
        filename = '%s.csv' % elem
        if os.stat(filename).st_size == 0:  # if file empty, write header
            print "%s is empty, adding header" % filename
            csv_dict[elem].writerow(header)
        else:
            print "%s already exists" % filename
    # append all data to csv
    for i, score_filename in enumerate(os.listdir(score_file_directory)):
        score_info = get_score_info(score_file_directory, score_filename)
        category = score_info[3]
        csv_dict[category].writerow(score_info)
    # close csvs
    for elem in categories:
        with open('%s.csv' % elem) as f:
            f.close()
    print "csv writing complete"


if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(description='Compile score files into a csv.')
    parser.add_argument('results_dir', metavar='D', type=str, help='directory of score files to collect')
    args = parser.parse_args()
    # write to csv
    results_to_csv(args.results_dir)
