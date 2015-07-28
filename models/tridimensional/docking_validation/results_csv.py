import argparse
import csv
import os

from constants import SCOREFILE_LINES, CSV_HEADER


categories = ['Chimera', '3DNA']  # tools used to create 64 or 256 pam variants


def get_pam_and_tool_from_filename(score_filename):
    """Extracts pam/tool info from a score filename
    Args:
        score_filename: filename
        score_file: specific score file name (including .txt)
    Returns:
        pam, pam_tool
    Example:
        Example score file: 'aggt_Chimera.txt'
        output: ['aggt', 'Chimera']
    """
    pam, pam_tool_txt = score_filename.split('_')
    pam_tool = pam_tool_txt.split('.')[0]
    return pam, pam_tool


def get_score_info(score_dir, score_filename):
    """Extracts pam/tool info and scores as a list from a score file
    Args:
        score_dir: path to score file directory
        score_file: specific score file name (including .txt)
    Returns:
        list of scores (as strings)
    Example:
        Example score file: 'aggt_Chimera.txt'
            Init FA Score:  2465.219
            Final FA Score:  2044.202 
            Init DNA score:  305.504
            Final DNA score:  -207.921
            Total time:  857.995
            Dock time:  600.995
        Expected return: ['a', 'g', 'g', 't', 'Chimera', '2465.219', ..., '600.995']
    """
    # read contents of textfile
    with open(os.path.join(score_dir, score_filename)) as f:
        file_lines = f.readlines()
        assert len(file_lines) == SCOREFILE_LINES
    # format scores
    return [file_lines[i].split(" ")[-1][:-1] for i in xrange(SCOREFILE_LINES)]


def results_to_csv(score_file_directory):
    """Appends data from each result file to an appropriate csv within score_file_directory
    Args:
        score_file_directory: path to score files
    Returns:
        None
    Notes:
        - currently creates / looks for a csv for each PAM tool
    """
    # csv prep
    csv_dict = {elem: csv.writer(open(os.path.join(score_file_directory, '%s.csv' % elem), 'a'), lineterminator='\n')
                for elem in categories}
    for elem in categories:
        filename = os.path.join(score_file_directory, '%s.csv' % elem)
        if os.stat(filename).st_size == 0:  # if file empty, write header
            print "%s is empty, adding header" % filename
            csv_dict[elem].writerow(CSV_HEADER)
        else:
            print "%s already exists" % filename
    # append all data to csv
    for i, score_filename in enumerate(os.listdir(score_file_directory)):
        if score_filename[-4:] == '.txt':
            pam, pam_tool = get_pam_and_tool_from_filename(score_filename)
            assert pam_tool in categories
            score_info = get_score_info(score_file_directory, score_filename)
            csv_dict[pam_tool].writerow(pam + pam_tool + score_info)
    # close csvs
    for elem in categories:
        with open(os.path.join(score_file_directory, '%s.csv' % elem), 'a') as f:
            f.close()
    print "csv writing complete"
    return


if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(description='Compile score files into a csv.')
    parser.add_argument('results_dir', metavar='D', type=str, help='directory of score files to collect')
    args = parser.parse_args()
    # write to csv
    results_to_csv(args.results_dir)
