# takes a genome FASTA and outputs a string of that genome

def convert(file_name):
	with open(file_name) as f:
		# splits by lines, removes first line, joins lines together
		sequence = ''.join(f.read().split('\n')[1:]).lower()
	return sequence