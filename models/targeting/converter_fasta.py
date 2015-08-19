def convert_fasta_to_string(filename):
    """Takes a genome FASTA and outputs a string of that genome
    Args:
        filename: fasta file
    Returns:
        string of the genome sequence
    """
    assert filename.split('.')[-1] == 'fasta'  # assert correct file type
    with open(filename) as f:
        sequence = ''.join(f.read().split('\n')[1:]).lower()  # splits by lines, removes first line, joins lines
    return sequence
