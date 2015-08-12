from random import randint

# random insertion function
# takes insertion size and produces random insertion string

DNA_ALPHABET = "acgt"

def nt_rand(insertion_size):
	insertion = ""
	for x in range(insertion_size):
		insertion += DNA_ALPHABET[randint(0,3)]

	return insertion
