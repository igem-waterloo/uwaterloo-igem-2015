import csv
import os


def results_to_csv(output_directory, csv_name, csv_header, data):
    """Write data to csv
    Args:
        data: a list of lists. Each element is a row
    Returns:
        None
    NOTE - suppressing output unconditionally for speedhax
    """
    # csv prep
    path_csv = os.path.join(output_directory, csv_name)
    writer_csv = csv.writer(open(path_csv, 'a'), lineterminator='\n')

    # if file empty, write header
    if os.stat(path_csv).st_size == 0:
        #print "%s is empty, adding header" % path_csv
        writer_csv.writerow(csv_header)
    #else:
    #    print "%s already exists" % path_csv

    # append all data to csv
    for csv_row in data:
        writer_csv.writerow(csv_row)
    # close csv
    with open(path_csv, 'a') as f:
        f.close()
    #print "csv writing complete"
    return


def csv_load(fullpath):
    """Read in csv and return dictionary of parsed data
    Args:
        fullpath: path to csv to be read
    Returns:
        parsed_dict: dictionary containing parsed data
    """
    assert fullpath[-4:] == '.csv'
    with open(fullpath, 'rb') as f:
        reader = csv.reader(f)
        csv_data = []
        for i, row in enumerate(reader):
            if i == 0:
                csv_header = row
            else:
                csv_data.append(row)
    return csv_header, csv_data


def csv_to_dict(fullpath):
    """Loads csv data into memory, then converts the data into dictionary format
    Args:
        fullpath: full path to file (e.g. root/some_dir/some_file.csv)
    Returns:
        csv in dictionary format where stats reference dictionary: {statname: {label: val_1 ... val_n} }
    Example dictionary:
        {'gene1': ["active",...,"deactivated"]
             ...
        }
    """
    csv_header, csv_data = csv_load(fullpath)
    column_index_dict = {key: csv_header.index(key) for key in csv_header}  # select columns for referencing data
    csv_dict = {}
    for key in csv_header:
        csv_dict[key] = [0]*len(csv_data)  # initialize list
        for i, row in enumerate(csv_data):
            csv_dict[key][i] = row[column_index_dict[key]]
    csv_dict['header'] = csv_header  # include header to preserve order
    return csv_dict
