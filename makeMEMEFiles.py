def makeMEMEFiles (PWMFileName, outputFileNamePrefix):
	# Convert PWM file into MEME file
	PWMFile = open(PWMFileName)
	PWMLines = PWMFile.readlines()
	PWMFile.close()
	TFLine = PWMLines[1]
	TFLineElements = TFLine.split("\t")
	TFName = TFLineElements[1].strip().upper()
	numPos = len(PWMLines) - 7
	PWMFileNameParts = PWMFileName.split("/")
	PWMFileNameElements = PWMFileNameParts[len(PWMFileNameParts)-1].split("_")
	outputFileName = outputFileNamePrefix + "_" + PWMFileNameElements[1] + "_" + PWMFileNameElements[2] + "_" + PWMFileNameElements[3]
	outputFile = open(outputFileName, 'w+')
	outputFile.write("MEME version 4.9.0\n\n")
	outputFile.write("ALPHABET= ACGT\n\n")
	outputFile.write("strands: +-\n\n")
	outputFile.write("MOTIF " + TFName + "\n")
	outputFile.write("letter-probability matrix: alength= 4 w= " + str(numPos) + "\n")
	for line in PWMLines[7:]:
		# Iterate through base probabilities at each position and write the probabilities to the output file
		lineElements = line.split("\t")
		for prob in lineElements[1:]:
			# Itereate through probabilities and write each to the outputFile
			outputFile.write(prob.strip())
			outputFile.write("\t")
		outputFile.write("\n")
	outputFile.close()

def makeMEMEFilesLoop (PWMFileNameListFileName, outputFileNamePrefix):
	# Convert all of the PWM files in a list into separate MEME files
	PWMFileNameListFile = open(PWMFileNameListFileName)
	for line in PWMFileNameListFile:
		# Iterate through the PWM files and create a MEME file for each
		makeMEMEFiles(line.strip(), outputFileNamePrefix)
	PWMFileNameListFile.close()

if __name__=="__main__":
   import sys
   PWMFileNameListFileName = sys.argv[1] 
   outputFileNamePrefix = sys.argv[2]
   makeMEMEFilesLoop (PWMFileNameListFileName, outputFileNamePrefix)
	
