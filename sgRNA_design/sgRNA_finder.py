from collections import defaultdict

# Function to determine if a potential PAM site meets the specified range and 
# efficiency values. pamLine should be one line of csv with format:
# Start Pos,Stop Pos, sgRNA sequence, PAM sequence (NGG), strand, MIT score
def pam_In_Range(pamLine,geneRange,effCutOff=0):

	# Parse and split csv line into list of string
	words = pamLine.strip("\n").strip(",").split(',')
	# Calculate the end position for PAM sequence based on desired coverage
	endRange = (geneRange[1] - geneRange[0])*geneRange[2] + geneRange[0]

	# Can be called directly in if and return statements, but makes code
	# difficult to read
	pamSta = int(words[0]) # Start Position
	pamStp = int(words[1]) # Stop Position
	pamSeq = words[2]      # sgRNA sequence
	pamNGG = words[3]      # PAM sequence
	pamStr = int(words[4]) # plus or minus strand
	pamEff = float(words[5]) # MIT efficiency score

	# determine if the PAM falls within the desired start and end sequence, if
	# it meets the desired efficiency cut off (set to 0 by default) and whether
	# it has masked characters in the sequence. Returns any PAMs that pass as 
	# a list of values. Returns empty string if it fails because None value
	# causes error in generating dictionary from returned values
	if pamSta > geneRange[0] and pamStp < endRange and 'N' not in pamSeq and \
	pamEff > effCutOff:
		return [pamSta,pamStp,pamSeq,pamNGG,pamStr,pamEff];
	else: return "";

# function to return dictionary with protein name as key and start & end
# positions as values. Takes csv file with format:
# Protein Name,Start Pos,End Pos,Coverage
def protein_Sites(fn):
	
	with open(fn, "r") as f:
		proteins = {}
		for line in f:
			# csv line formatting
			words = line.strip("\n").strip(",").split(',')
			proteins[words[0]] = (int(words[1]),int(words[2]),float(words[3]))
	f.close()
	return proteins;

# Sort list of PAMs (which are list themselves) into ordered list based on start
# position. Starting with first PAM, determines if two sites overlap (with 
# extension for Cas9 size) and discards the lower scoring one. Is greedy, and 
# will discard two non-overlapping PAMs for a third that overlaps both if that
# one is higher scoring.
def pam_Overlap(listOfPams,overlapExtension=0):

	# I have no idea how this works, but it gives a sorted list of PAMs. Thank
	# you StackOverflow
	# stackoverflow.com/questions/5201191/sorting-a-list-of-lists-in-python
	listOfPams.sort(key=lambda x: x[1])	
	i=0
	while i+1 < len(listOfPams):
		firstPAM = listOfPams[i]
		secondPAM = listOfPams[i+1]

		if secondPAM[0] - firstPAM[0] < 23 + overlapExtension:
			if firstPAM[5] < secondPAM[5]:
				listOfPams.pop(i)
				#print("Discarded:" + str(firstPAM))
			else: 
				#print("Discarded:" + str(secondPAM))
				listOfPams.pop(i+1)

		else:
			i+=1

if __name__ == '__main__':
	f1 = input('Enter name of gene range file:\n')
	
	proteins = protein_Sites(f1)
	print(proteins)

	f2 = input('Enter name of possible PAM site file:\n')

	pamListByProtein = defaultdict(list)

	with open(f2, "r") as f:
		for line in f:
			for key in proteins:
				pam = pam_In_Range(line,proteins[key])
				if len(pam) > 1: pamListByProtein[key].append(pam)

	f.close()

	with open("output.csv", "w") as output:

		for key in sorted(pamListByProtein):
			pam_Overlap(pamListByProtein[key])
			for i,v in enumerate(pamListByProtein[key]):
				print(str(v[0])+',' + str(v[1])+',' + v[2]+',' + v[3]+',' + str(v[4])+',' + str(v[5]), file = output)

	output.close()
