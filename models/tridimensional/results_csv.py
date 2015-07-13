import argparse
import csv
import os


# helper function for extracting score data
def get_scores(filedata):
    # Given a score file, extract the last element of each row
    # (expected to be a float) and remove the last 2 chars (a newline)
    assert len(filedata) == 6
    return [filedata[i].split(" ")[-1][:-1] for i in xrange(4)]


# argument parsing
parser = argparse.ArgumentParser(description='Compile score files into a csv.')
parser.add_argument('results_dir', metavar='D', type=str, 
                   help='directory of score files to collect')
args = parser.parse_args()
score_file_directory = args.results_dir
score_files = os.listdir(score_file_directory)  # get list of score files

# csv prep
header = ['PAM_1','PAM_2','PAM_3','Tool','Init FA','Init DNA','Final FA','Final DNA']
categories = ['3DNA', 'Chimera']
csv_dict = {elem:csv.writer(open('%s.csv' % elem,'a'),lineterminator='\n') for elem in categories}
for elem in categories:
    filename = '%s.csv' % elem
    if os.stat(filename).st_size == 0:  # if file empty, write header
        print "%s is empty, adding header" % filename
        csv_dict[elem].writerow(header)
    else:
        print "%s already exists" % filename
        
# append all data to csv
for i, score_file in enumerate(score_files):
    # get textfile classifiers (PAM and category)
    pam, pam_tool_txt = score_file.split("_")
    pam_tool = pam_tool_txt.split(".")[0]
    assert pam_tool in categories
    # read contents of textfile
    with open(score_file_directory + score_file) as f:
        filedata = f.readlines()
    # format row and append to csv
    scores = get_scores(filedata)
    row = [pam[0], pam[1], pam[2], pam_tool] + scores
    csv_dict[pam_tool].writerow(row)

# close csvs
for elem in categories:
    with open('%s.csv' % elem) as f: f.close()
print "csv writing complete"
