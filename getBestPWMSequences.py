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
	PWMFileNameList = open(PWMFileNameListFile)
	bestPWMDict = {}
	TFList = []
	for PWMFileName in PWMFileNameList:
		# Iterate through PWMs and find the best sequence for each
		PWMFileNamePathElements = PWMFileName.strip().split("/")
		PWMFileNameElements = PWMFileNamePathElements[len(PWMFileNamePathElements)-1].strip().split("_")
		TF = PWMFileNameElements[1].upper()
		if TF in bestPWMDict:
			# ALWAYS CHOOSING THE FIRST TF
			continue
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
		bestPWMDict[TF] = bestPWM
		TFList.append(TF)
	PWMFileNameList.close()
	return [bestPWMDict, TFList]


def compareBestPWMs(PWMFileNameListFileOne, PWMFileNameListFileTwo):
	# Compare the best PWMs for TFs between two different species
	[bestPWMDictOne, TFListOne] = getBestPWMSequences(PWMFileNameListFileOne)
	[bestPWMDictTwo, TFListTwo] = getBestPWMSequences(PWMFileNameListFileTwo)
	commonTFList = []
	commonFracList = []
	for TF in TFListOne:
		# Iterate through TFs and determine the fraction of bases in the PWMs that are the same for each TF across species
		if TF in TFListTwo:
			# The TF is in both species
			bestPWMOne = bestPWMDictOne[TF]
			bestPWMTwo = bestPWMDictTwo[TF]
			smaller = bestPWMOne
			larger = bestPWMTwo
			if len(bestPWMOne) > len(bestPWMTwo):
				# The second PWM is shorter
				smaller = bestPWMTwo
				larger = bestPWMOne
			if len(smaller) == 0:
				# The smaller TF has an empty motif, so continue
				continue
			commonTFList.append(TF)
			commonNumBasesList = []
			for i in range(len(larger)-len(smaller)+1):
				# Iterate through positions in the larger PWM where the smaller PWM might align
				commonNumBases = 0
				for j in range(len(smaller)):
					# Iterate through bases and count the number that are the same in both PWMs
					if larger[i+j] == smaller[j]:
						# The bases are the same at the current location in the best PWMs
						commonNumBases = commonNumBases + 1
				commonNumBasesList.append(commonNumBases)
			commonFrac = float(max(commonNumBasesList))/float(len(smaller))
			commonFracList.append(commonFrac)
	return [commonTFList, commonFracList]


def writeTFInfoToFile(commonTFList, commonFracList, outputFileName):
	# Write the fraction of common bases across species for each TF to a file
	outputFile = open(outputFileName, 'w+')
	for i in range(len(commonTFList)):
		# Iterate through common TFs across the species and write each TF's information to the output file
		outputFile.write(commonTFList[i])
		outputFile.write("\t")
		outputFile.write(str(commonFracList[i]))
		outputFile.write("\n")
	outputFile.close()


if __name__=="__main__":
   import sys
   PWMFileNameListFileOne = sys.argv[1] 
   PWMFileNameListFileTwo = sys.argv[2]
   outputFileName = sys.argv[3]

   [commonTFList, commonFracList] = compareBestPWMs(PWMFileNameListFileOne, PWMFileNameListFileTwo)
   writeTFInfoToFile(commonTFList, commonFracList, outputFileName)
