import argparse
import csv
import os
import shutil
from tempfile import NamedTemporaryFile

from constants import CSV_HEADER, CSV_HEADER_64, DNA_ALPHABET, PAM_TOOLS, SCOREFILE_LINES


def get_pam_and_tool_from_filename(score_filename):
    """Extracts pam/tool info from a score filename
    Args:
        score_filename: specific score file name (including .txt)
    Returns:
        list of pam, pam_tool
    Example:
        Example score file: 'aggt_Chimera.txt'
        output: ['aggt', 'Chimera']
    """
    pam, pam_tool_txt = score_filename.split('_')
    pam_tool = pam_tool_txt.split('.')[0]
    return pam, pam_tool


def get_score_info(score_dir, score_filename):
    """Extracts pam/tool info and scores as a list from a score fileW
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
        Expected return: ['2465.219', ..., '600.995']
    """
    # read contents of textfile
    with open(os.path.join(score_dir, score_filename)) as f:
        file_lines = f.readlines()
        assert len(file_lines) == SCOREFILE_LINES
    # format scores
    return [file_lines[i].split(" ")[-1][:-1] for i in xrange(SCOREFILE_LINES)]


def csv_header_bugfix(fullpath):
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
        print "Warning: this csv has a bug -- columns expect 4 letter pam, data is for 3 letter pam"
        print "Notice: fixing the file by removing 'PAM_4' from header"
        tempfile = NamedTemporaryFile(delete=False)
        with open(fullpath, 'rb') as csvfile, tempfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            writer = csv.writer(tempfile, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    row.pop(3)
                writer.writerow(row)
        shutil.move(tempfile.name, fullpath)
    return


def results_to_csv(score_file_directory, pam_tool='Chimera'):
    """Appends data from each result file to an appropriate csv within score_file_directory
    Args:
        score_file_directory: path to score files
        pam_tool: [string] either "3DNA" or "Chimera"
    Returns:
        None
    """
    # csv prep
    assert pam_tool in PAM_TOOLS
    path_csv = os.path.join(score_file_directory, '%s.csv' % pam_tool)
    score_writer = csv.writer(open(path_csv, 'a'), lineterminator='\n')
    # get pam length for header condition
    score_filelist = os.listdir(score_file_directory)
    for score_filename in score_filelist:
        if score_filename[-4:] == '.txt' and pam_tool in score_filename:
            pam = get_pam_and_tool_from_filename(score_filename)[0]
            if len(pam) == 3:
                csv_header = CSV_HEADER_64
            else:
                csv_header = CSV_HEADER
            break
    # if file empty, write header
    if os.stat(path_csv).st_size == 0:
        print "%s is empty, adding header" % path_csv
        score_writer.writerow(csv_header)
    else:
        print "%s already exists" % path_csv
    # append all data to csv
    for i, score_filename in enumerate(score_filelist):
        if score_filename[-4:] == '.txt' and pam_tool in score_filename:
            pam, tool = get_pam_and_tool_from_filename(score_filename)
            assert len(pam) == 3 or len(pam) == 4
            assert tool == pam_tool
            score_info = get_score_info(score_file_directory, score_filename)
            csv_row = list(pam) + [pam_tool] + score_info  # all elements must be lists to append them
            score_writer.writerow(csv_row)
    # close csv
    with open(path_csv, 'a') as f:
        f.close()
    print "csv writing complete"
    return


if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(description='Compile score files into a csv.')
    parser.add_argument('--results_dir', metavar='D', type=str, help='directory of score files to collect')
    parser.add_argument('--alt_tool', metavar='S', nargs='?', const='3DNA', default='Chimera',
                        type=str, help='[switch] compile 3DNA score files (default: Chimera)')
    args = parser.parse_args()
    # write to csv
    results_to_csv(args.results_dir, pam_tool=args.alt_tool)
