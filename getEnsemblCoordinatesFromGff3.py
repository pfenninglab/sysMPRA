import sys
import argparse

def parseArgument():
	# Parse the input
	parser=argparse.ArgumentParser(description=\
			"Get Ensembl gene coordinates from a gff3 file")
	parser.add_argument("--gff3FileName", required=True,\
		help='gff3 file')
	parser.add_argument("--genesFileName", required=False, default=None,\
		help='File with list of genes, if None then use all of the genes')
	parser.add_argument("--geneName", action="store_true", required=False,\
		help='Use gene names instead of Ensembl gene names')
	parser.add_argument("--geneNameIndex", type=int, required=False, default=5,\
                help='Index in description of gene name; 0-indexed')
	parser.add_argument("--convertToUpper", action="store_true", \
		required=False,\
                help='Convert each gene name to upper-case')
	parser.add_argument("--outputFileName", required=True,\
		help='File where gene coordinates will be recorded')
	options = parser.parse_args()
	return options

def getEnsemblCoordinatesFromGff3(options):
	# Get Ensembl gene coordinates from a gff3 file
	gff3File = open(options.gff3FileName)
	genesFile = None
	genes = None
	if options.genesFileName != None:
		# A file with a subset of the genes has been provided
		genesFile = open(options.genesFileName)
		genes = [line.strip() for line in genesFile]
		if options.convertToUpper:
			# Convert the gene names to upper-case
			genes = [g.upper() for g in genes]
		genesFile.close()
	outputFile = open(options.outputFileName, 'w+')
	for line in gff3File:
		# Iterate through the lines of the gff3 file and save the coordinates of the genes in the list
		lineElements = line.strip().split("\t")
		attributeElements = lineElements[8].split(";")
		if options.geneName:
			# Use the gene names instead of Ensembl gene names
			geneNameElements =\
				attributeElements[options.geneNameIndex].split("=")
			currentGeneName = geneNameElements[1]
			if options.convertToUpper:
				# Convert the gene name to upper-case
				currentGeneName = currentGeneName.upper()
			if (not options.genesFileName) or\
				(currentGeneName in genes):
				# The current gene is in the list, so record its information
				outputFile.write("\t".join([lineElements[0], lineElements[3], \
					lineElements[4], currentGeneName, \
					lineElements[6], lineElements[5]]) + \
					"\n")
		else:
			# Use Ensembl gene names
			geneIdElements = attributeElements[2].split("=")
			geneId = geneIdElements[1].split(".")[0]
			if (not options.genesFileName) or (geneId in genes):
				# The current gene is in the list, so record its information
				outputFile.write("\t".join([lineElements[0], \
					lineElements[3], lineElements[4], \
					geneId, lineElements[6], \
					lineElements[5]]) + "\n")
	gff3File.close()
	outputFile.close()

if __name__ == "__main__":
	options = parseArgument()
	getEnsemblCoordinatesFromGff3(options)
