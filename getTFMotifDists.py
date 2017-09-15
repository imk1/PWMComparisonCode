def makeTFMotifFileListArrayPair(TFMotifFileListFileName):
	# Iterate through the TF motif files and make a pair of arrays that maps TF motifs to TF motif file names
	TFMotifFileListFile = open(TFMotifFileListFileName)
	TFMotifFileListTFMotifs = []
	TFMotifFileListFileNames = []
	for line in TFMotifFileListFile:
		# Iterate through the file names and enter each in the entry of the appropriate TF motif
		lineElements = line.strip().split("/")
		nameElements = lineElements[-2].split("_")
		TF = nameElements[1].upper()
		motif = nameElements[2]
		if (TF, motif) in TFMotifFileListTFMotifs:
			# There is already an entry for the TF, so append the file name to that entry
			TFMotifIndex = TFMotifFileListTFMotifs.index((TF, motif))
			TFMotifFileListFileNames[TFMotifIndex].append(line.strip())
		else:
			TFMotifFileListTFMotifs.append((TF, motif))
			TFMotifFileListFileNames.append([line.strip()])
	TFMotifFileListFile.close()
	return [TFMotifFileListTFMotifs, TFMotifFileListFileNames]

def getTFMotifDists(TomtomFileListFileName, outputFileName):
	# Get the distance from each motif to the closest motif in the other species for the same TF
	# HIGHER DISTANCE MEANS CLOSER
	# Output will be TF, motif from species 1, best motif from species 2 for motif from species 1, motif distance
	[TomtomFileListTFMotifs, TomtomFileListFileNames] = makeTFMotifFileListArrayPair(TomtomFileListFileName)
	outputFile = open(outputFileName, 'w+')
	for i in range(len(TomtomFileListTFMotifs)):
		# Iterate through the Tomtom files and write the best for each motif pair from the first species
		bestDist = 0 - float('inf')
		bestMotif = ""
		for TomtomFileName in TomtomFileListFileNames[i]:
			# Iterate through the file names for the motifs mapped to the (TF, motif) from the other species
			TomtomFile = open(TomtomFileName)
			TomtomFileLines = TomtomFile.readlines()
			TomtomFile.close()
			if len(TomtomFileLines) < 2:
				# The Tomtom file is blank, so a motif must be empty
				continue
			TomtomDistElements = TomtomFileLines[1].split("\t")
			TomtomDist = float(TomtomDistElements[3])
			TomtomOverlap = float(TomtomDistElements[7])
			TomtomDistNormalized = TomtomDist/TomtomOverlap # NORMALIZING SCORES BASED ON OVERLAP LENGTH
			if TomtomDistNormalized > bestDist:
				# The current distance is the best so far
				bestDist = TomtomDistNormalized
				lineElements = TomtomFileName.strip().split("/")
				nameElements = lineElements[-2].split("_")
				motifTwo = nameElements[3]
				bestMotif = motifTwo
		if bestDist == 0 - float('inf'):
			# The motif in the first species or all motifs in the second species are empty, so continue
			continue
		outputFile.write(TomtomFileListTFMotifs[i][0] + "\t" + TomtomFileListTFMotifs[i][1] + "\t" + bestMotif + "\t" + str(bestDist) + "\n")
	outputFile.close()

if __name__=="__main__":
   import sys
   TomtomFileListFileName = sys.argv[1] 
   outputFileName = sys.argv[2]
   getTFMotifDists(TomtomFileListFileName, outputFileName)
