def make_pwm_cisbp(txt,pseudo_count):
   #takes an input PWM file (in cis-bp format) and outputs a PWM array ...
   #txt: is the name of the PWM file (7 descriptor rows, 1 column with numbers and 4 columns for ACGT respectively, the number additional rows corresponds ...
   #to the length of the motfif
   #psuedo_count can be set to a small value (e.g. 0.001) to avoid over-fitting       

   o = open(txt)
   lines = o.readlines()[7:] # Remove descriptor lines
   o.close()

   seq = 'ACGT'
   mat = {}
   total_scores = {}
    
   for j in range(0,len(lines)):
       # Iterate through positions and compute the total score for each position
       total_scores[j] = 0
       for i in range(0,4):
	   # Add the score for each base to the total score
           total_scores[j] +=  float(lines[j].strip().split()[i+1])

   for i in range(0,4):
      # Compute the score for each base at each position by dividing the score for the base by the total score and incorporating pseudocounts
      mat[i] = []     
      for j in range(0,len(lines)):
	 # Iterate through positions and compute the score for each position
	 s = lines[j].strip().split()
         mat[i].append((float(pseudo_count + float(s[i+1])))/float(4*pseudo_count + total_scores[j]))
   return mat


def getBestPWMSequences(PWMFileNameListFile):
	# Creates a dictionary with the sequences for the best PWMs for all TFs
	# ALLOWS FOR MULTIPLE PWMS FOR EACH TF
	PWMFileNameList = open(PWMFileNameListFile)
	bestPWMDict = {}
	TFList = []
	for PWMFileName in PWMFileNameList:
		# Iterate through PWMs and find the best sequence for each
		PWMFileNamePathElements = PWMFileName.strip().split("/")
		PWMFileNameElements = PWMFileNamePathElements[len(PWMFileNamePathElements)-1].strip().split("_")
		TF = PWMFileNameElements[1].upper()
		mat = make_pwm_cisbp(PWMFileName.strip(), .0001)
		bestPWM = []
		for i in range(len(mat[0])):
			# Iterate through positions and find the best base for each position
			currentBestNum = 0
			currentBestProb = mat[0][i]
			for j in range(1, 4):
				# Iterate through bases to find the most likely base
				if mat[j][i] > currentBestProb:
					# A more likely base has been identified, so set the best to be this base
					# IF 2 BASES HAVE THE SAME PROBABILITY, ALWAYS CHOOSES THE FIRST BASE
					currentBestNum = j
					currentBestProb = mat[j][i]
			bestPWM.append(currentBestNum)
		if len(bestPWM) == 0:
			# The PWM file is empty, so continue
			continue
		if TF not in bestPWMDict:
			# Initialize the dictionary for the current TF
			bestPWMDict[TF] = []
			TFList.append(TF)
		bestPWMDict[TF].append(bestPWM)
	PWMFileNameList.close()
	return [bestPWMDict, TFList]


def compareBestPWMs(PWMFileNameListFileOne, PWMFileNameListFileTwo):
	# Compare the best PWMs for TFs between two different species
	# If there are multiple PWMs for the current TF, find the closest one for each
	[bestPWMDictOne, TFListOne] = getBestPWMSequences(PWMFileNameListFileOne)
	[bestPWMDictTwo, TFListTwo] = getBestPWMSequences(PWMFileNameListFileTwo)
	commonTFList = []
	commonFracList = []
	for TF in TFListOne:
		# Iterate through TFs and determine the fraction of bases in the PWMs that are the same for each TF across species
		if TF in TFListTwo:
			# The TF is in both species
			bestPWMOneList = bestPWMDictOne[TF]
			bestPWMTwoList = bestPWMDictTwo[TF]
			for bestPWMOne in bestPWMOneList:
				# Iterate through the best motifs for the first species
				currentBestFrac = 0
				for bestPWMTwo in bestPWMDictTwo[TF]:
					# Iterate through the best motifs in the second species for the current TF
					smaller = bestPWMOne
					larger = bestPWMTwo
					if len(bestPWMOne) > len(bestPWMTwo):
						# The second PWM is shorter
						smaller = bestPWMTwo
						larger = bestPWMOne
					commonNumBasesList = []		
					for i in range(len(larger)-len(smaller)+1):
						# Iterate through positions in larger PWM where smaller PWM may align
						commonNumBases = 0
						for j in range(len(smaller)):
							# Iterate through bases and count number that are same in both PWMs
							if larger[i+j] == smaller[j]:
								# Bases are the same at current location in best PWMs
								commonNumBases = commonNumBases + 1
						commonNumBasesList.append(commonNumBases)
					commonFrac = float(max(commonNumBasesList))/float(len(smaller))
					if commonFrac > currentBestFrac:
						# The current motif aligns better than the best motif so far
						currentBestFrac = commonFrac
				commonFracList.append(currentBestFrac)
				commonTFList.append(TF)
	return [commonTFList, commonFracList]

def makeTFFamilyDict(TFInformationFileName):
	# Make a dictionary that maps TFs to families
	TFFamilyDict = {}
	TFInformationFile = open(TFInformationFileName)
	TFInformationFile.readline() # Remove the header
	for line in TFInformationFile:
		# Iterate through the information about each TF and enter each TF and its family into the dictionary
		lineElements = line.split("\t")
		TFFamilyDict[lineElements[1].upper()] = lineElements[5]
	TFInformationFile.close()
	return TFFamilyDict


def writeTFInfoToFile(commonTFList, commonFracList, TFFamilyDict, outputFileName):
	# Write the fraction of common bases across species for each TF to a file
	outputFile = open(outputFileName, 'w+')
	for i in range(len(commonTFList)):
		# Iterate through common TFs across the species and write each TF's information (name, family, similarity fraction for current motif) to the output file
		outputFile.write(commonTFList[i])
		outputFile.write("\t")
		outputFile.write(TFFamilyDict[commonTFList[i]])
		outputFile.write("\t")
		outputFile.write(str(commonFracList[i]))
		outputFile.write("\n")
	outputFile.close()


if __name__=="__main__":
	import sys
	PWMFileNameListFileOne = sys.argv[1] 
	PWMFileNameListFileTwo = sys.argv[2]
	TFInformationFileName = sys.argv[3]
	outputFileName = sys.argv[4]
	[commonTFList, commonFracList] = compareBestPWMs(PWMFileNameListFileOne, PWMFileNameListFileTwo)
	TFFamilyDict = makeTFFamilyDict(TFInformationFileName)
	writeTFInfoToFile(commonTFList, commonFracList, TFFamilyDict, outputFileName)
