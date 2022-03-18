import sys
import argparse
import numpy as np

def parseArgument():
	# Parse the input
	parser=argparse.ArgumentParser(description=\
			"Find the best scrambled motif, meaning the motif with the weakest other motif hit")
	parser.add_argument("--fimoBedFileName", required=True,\
			help='Name of the bed file with fimo results, sorted by region, then scramble, then motif hit strength')
	parser.add_argument("--outputFileName", required=True,\
			help='File with the best scramble')
	parser.add_argument("--numScrambles", type=int, required=False, default=500,\
			help='File with the best scramble')
	options = parser.parse_args()
	return options

def chooseBestScramble(options):
	# Find the best scrambled motif, meaning the motif with the weakest other motif hit
	fimoFile = open(options.fimoBedFileName)
	outputFile = open(options.outputFileName, 'w+')
	lastRegion = ""
	lastScramble = ""
	bestScramble = ""
	bestpVal = 0.0
	scrambleHasMotif = np.zeros(options.numScrambles)
	for line in fimoFile:
		# Iterate through the scramble, p-value combinations and get the best scramble (the scramble with the highest lowest p-value) for each
		lineElements = line.strip().split("\t")
		scrambleElements = lineElements[3].split("_")
		currentRegion = scrambleElements[0]
		scrambleNum = int(scrambleElements[2]) - 1
		if lastRegion == "":
			# At the beginning
			lastRegion = currentRegion
		if currentRegion != lastRegion:
			# At a new region, so record the best scramble and re-set everything
			if np.sum(scrambleHasMotif) < options.numScrambles:
				# There is a scramble with no motif hit
				bestScrambleIndex = np.argmin(scrambleHasMotif)[0]
				bestScramble = "_".join([currentRegion, "shuf", str(bestScrambleIndex)])
				bestpVal = 1.0
			outputFile.write("\t".join([lastRegion, bestScramble, str(bestpVal)]) + "\n")
			lastScramble = ""
			bestScramble = ""
			bestpVal = 0.0
			lastRegion = currentRegion
			scrambleHasMotif = np.zeros(options.numScrambles) # Accidentally excluded when designing Mef2c scrambles for MPRA, but did not affect results
		scrambleHasMotif[scrambleNum] = 1
		if lastScramble != lineElements[3]:
			# At a new scramble, so check if the new lowest p-value is higher than the previous oldest p-value
			currentpVal = float(lineElements[4])
			if currentpVal > bestpVal:
				# The p-value for the current scramble is the best one
				bestpVal = currentpVal
				bestScramble = lineElements[3]
			lastScramble = lineElements[3]
	outputFile.write("\t".join([lastRegion, bestScramble, str(bestpVal)]) + "\n")
	fimoFile.close()
	outputFile.close()

if __name__ == "__main__":
	options = parseArgument()
	chooseBestScramble(options)
