import csv
import os


def results_to_csv(output_directory, csv_name, csv_header, data):
    """Write data to csv
    Args:
        data: a list of lists. Each element is a row
    Returns:
        None
    """
    # csv prep
    path_csv = os.path.join(output_directory, csv_name)
    writer_csv = csv.writer(open(path_csv, 'a'), lineterminator='\n')

    # if file empty, write header
    if os.stat(path_csv).st_size == 0:
        print "%s is empty, adding header" % path_csv
        writer_csv.writerow(csv_header)
    else:
        print "%s already exists" % path_csv

    # append all data to csv
    for csv_row in data:
        writer_csv.writerow(csv_row)
    # close csv
    with open(path_csv, 'a') as f:
        f.close()
    print "csv writing complete"
    return
