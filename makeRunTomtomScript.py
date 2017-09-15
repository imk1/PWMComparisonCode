def makeTFFileListArrayPair(TFFileListFileName):
	# Iterate through the TF files and make a pair of arrays that maps TFs to TF file names
	TFFileListFile = open(TFFileListFileName)
	TFFileListTFs = []
	TFFileListFileNames = []
	for line in TFFileListFile:
		# Iterate through the file names and enter each in the entry of the appropriate TF
		lineElements = line.strip().split("/")
		nameElements = lineElements[-1].split("_")
		TF = nameElements[1].upper()
		if TF in TFFileListTFs:
			# There is already an entry for the TF, so append the file name to that entry
			TFIndex = TFFileListTFs.index(TF)
			TFFileListFileNames[TFIndex].append(line.strip())
		else:
			TFFileListTFs.append(TF)
			TFFileListFileNames.append([line.strip()])
	TFFileListFile.close()
	return [TFFileListTFs, TFFileListFileNames]

def makerunTomtomScript(MEMEFileListOneFileName, MEMEFileListTwoFileName, outputFileNamePrefix, dist, scriptFileName, queryPseudo, targetPseudo, icCutoff):
	# Makes a script that will run Tomtom for all pairs of motifs from the same TF
	MEMEFileListOneFile = open(MEMEFileListOneFileName)
	[MEMEFileListTwoTFs, MEMEFileListTwoFileNames] = makeTFFileListArrayPair(MEMEFileListTwoFileName)
	scriptFile = open(scriptFileName, 'w+')
	for line in MEMEFileListOneFile:
		# Iterate through the MEME files and write a command for Tomtom for each TF pair from the different species
		lineElements = line.strip().split("/")
		nameElements = lineElements[-1].split("_")
		TF = nameElements[1].upper()
		if TF in MEMEFileListTwoTFs:
			# The TF is present in both species
			TFIndex = MEMEFileListTwoTFs.index(TF)
			motifOne = nameElements[2]
			for fileName in MEMEFileListTwoFileNames[TFIndex]:
				# Iterate through the file names for the motifs with the TF from the other species
				lineElementsTwo = fileName.strip().split("/")
				nameElementsTwo = lineElementsTwo[-1].split("_")
				motifTwo = nameElementsTwo[2]
				outputFileName = outputFileNamePrefix + "_" + TF + "_" + motifOne + "_" + motifTwo
				inputFileList = line.strip() + " " + fileName
				scriptFile.write("tomtom -o " + outputFileName + " -dist " + dist + " -incomplete-scores -internal -query-pseudo " + queryPseudo + " -target-pseudo " + targetPseudo + " -ic-cutoff " + icCutoff + " " + inputFileList + "\n")
	MEMEFileListOneFile.close()
	scriptFile.close()

if __name__=="__main__":
	import sys
	MEMEFileListOneFileName = sys.argv[1] 
	MEMEFileListTwoFileName = sys.argv[2]
	outputFileNamePrefix = sys.argv[3]
	dist = sys.argv[4]
	scriptFileName = sys.argv[5]
	queryPseudo = sys.argv[6] # Should be a float
	targetPseudo = sys.argv[7] # Should be a float
	icCutoff = sys.argv[8] # Should be a float
	makerunTomtomScript(MEMEFileListOneFileName, MEMEFileListTwoFileName, outputFileNamePrefix, dist, scriptFileName, queryPseudo, targetPseudo, icCutoff)
