def convert_fasta_to_string(file_name):
    # takes a genome FASTA and outputs a string of that genome
    with open(file_name) as f:
        # splits by lines, removes first line, joins lines together
        sequence = ''.join(f.read().split('\n')[1:]).lower()
    return sequence
