def getPWMInformationContent(PWMMEMEFileName, pseudoCount):
	# Gets the information content of a PWM
	# Information content = self-information = -sum_positions(sum_bases(p_position,base * log_2(p_position,base)))
	PWMMEMEFile = open(PWMMEMEFileName)
	for i in range(8):
		# Assumes that the first 8 lines are header
		PWMMEMEFile.readline()
	selfInformation = 0
	for line in PWMMEMEFile:
		# Iterate through positions and add the self-information for each position to the total self-information
		lineElements = line.strip().split("\t")
		for pijstr in lineElements:
			# Iterate through the bases and add the self-information for each base to the total self-information
			pij = float(pijstr.strip()) + pseudoCount # NEED TO NORMALIZE BY TOTAL PROB. OF POSITION + 4*PSEUDOC
			selfInformationpij = 0 - (pij * numpy.log2(pij))
			selfInformation = selfInformation + selfInformationpij
	PWMMEMEFile.close()
	return selfInformation

def getPWMInformationContentLoop(PWMMEMEFileListFileName, pseudoCount, outputFileName):
	# Gets the information content for a list of PWMs
	# Each line of the output file is TF, motif, information content
	PWMMEMEFileListFile = open(PWMMEMEFileListFileName)
	outputFile = open(outputFileName, 'w+')
	for line in PWMMEMEFileListFile:
		# Iterate through the PWMs and get the information content for each
		selfInformation = getPWMInformationContent(line.strip(), pseudoCount)
		lineElements = line.split("/")
		nameElements = lineElements[-1].strip().split("_")
		outputFile.write(nameElements[1] + "\t" + nameElements[2] + "\t" + str(selfInformation) + "\n")
	PWMMEMEFileListFile.close()
	outputFile.close()

if __name__=="__main__":
   import sys
   import numpy
   PWMMEMEFileListFileName = sys.argv[1] 
   outputFileName = sys.argv[2]
   pseudoCount = float(sys.argv[3]) # Used to prevent 0 probabilities
   getPWMInformationContentLoop(PWMMEMEFileListFileName, pseudoCount, outputFileName)
