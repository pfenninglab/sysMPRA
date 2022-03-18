import sys
import argparse
from Bio import SeqIO

def parseArgument():
	# Parse the input
	parser=argparse.ArgumentParser(description=\
			"Filter a fasta file based on IDs")
	parser.add_argument("--fastaFileName", required=True,\
			help='Name of the original fasta file')
	parser.add_argument("--fastaIDsFileName", required=True,\
			help='Name of file with fasta IDs that will be selected')
	parser.add_argument("--IDCol", type=int, required=False, default=1, \
			help='Number of column with fasta IDs, zero-indexed')
	parser.add_argument("--outputFileName", required=True,\
			help='fasta file where the selected sequences will be written')
	options = parser.parse_args()
	return options

def filterFastaFile(options):
	# Filter a fasta file based on IDs
	fastaIDsFile = open(options.fastaIDsFileName)
	fastaIDs = []
	for line in fastaIDsFile:
		# Iterate through the fasta IDs and record a line for each
		lineElements = line.strip().split("\t")
		fastaIDs.append(lineElements[options.IDCol])
	fastaIDsFile.close()
	outputFile = open(options.outputFileName, 'w+')
	for seqRecord in SeqIO.parse(options.fastaFileName, "fasta"):
		# Iterate through the lines of the fasta file and select those with an ID on the list
		if seqRecord.id in fastaIDs:
			# Record the current fasta entry
			outputFile.write(">" + seqRecord.id + "\n")
			outputFile.write(str(seqRecord.seq) + "\n")
	outputFile.close()

if __name__ == "__main__":
	options = parseArgument()
	filterFastaFile(options)
