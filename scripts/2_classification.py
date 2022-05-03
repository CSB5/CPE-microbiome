import numpy as np
import os
from distutils.dir_util import copy_tree
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def listdir_nohidden(path):
    return [f for f in sorted(os.listdir(path)) if not f.startswith('.')]

#########################################################################################################
# Find coverage of the samples we have, do classification
#########################################################################################################

'''
# Change to kp as required
strain = "ec"
strain_path = "vcf_filtered_Jon_sorted/ec_filtered_vcf/"

tpList = listdir_nohidden(strain_path)	

# (5a) silhouette score for 2 class case
silScoreMatrix = np.zeros((sampleCount + 1, 5), dtype = object)
silScoreMatrix[0, 0:4] = ["indiv", "fam", "relation", "timepoint"]
silScoreMatrix[0, 4] = "2Class"

# (5b) class centers for 2 class case
cenClassMatrix = np.zeros((sampleCount + 1, 5), dtype = object)
cenClassMatrix[0, 0:4] = ["indiv", "fam", "relation", "timepoint"]
cenClassMatrix[0, 4] = "2Class"

# (5c) histogram bin counts (0 to 1 in 0.01 intervals)
histBinMatrix = np.zeros((sampleCount + 1, 104), dtype = object)
histBinMatrix[0, 0:4] = ["indiv", "fam", "relation", "timepoint"]
histBinMatrix[0, 4:104] = np.divide(range(100), 100.0)

print (histBinMatrix[0, :])

count = 1

for timepoint in tpList:

	vcfNoDup = np.loadtxt(strain_path + timepoint, dtype = str)

	print vcfNoDup

	cenClassMatrix[count, 0] = timepoint.split("_")[0]
	cenClassMatrix[count, 1:4] = timepoint.split("_")[0].split("-")

	silScoreMatrix[count, 0] = timepoint.split("_")[0]
	silScoreMatrix[count, 1:4] = timepoint.split("_")[0].split("-")

	histBinMatrix[count, 0] = timepoint.split("_")[0]
	histBinMatrix[count, 1:4] = timepoint.split("_")[0].split("-")

	afVec = vcfNoDup[:, 7].astype(float).reshape(-1, 1)

	# Filter for SNPs with <0.98 AF to do k-means clustering 
	segIndex = np.where(afVec < 0.98)[0]
	segAFVec = afVec[segIndex]

	km = KMeans(n_clusters = 2, random_state = 976).fit(segAFVec)
	preds = km.fit_predict(segAFVec)

	centStr = ""
	for cent in km.cluster_centers_:
		centStr = centStr + str(cent[0]) + ","

	cenClassMatrix[count, 4] = centStr
	silScoreMatrix[count, 4] = silhouette_score(segAFVec, preds, sample_size = 10000, random_state = 976)

	histBinMatrix[count, 4:104] = np.histogram(afVec, np.divide(range(101), 100.0))[0]

	count += 1

np.savetxt("5a_" + strain + "_sil_score.txt", silScoreMatrix, fmt = "%s", delimiter = "\t")
np.savetxt("5b_" + strain + "_class_cen.txt", cenClassMatrix, fmt = "%s", delimiter = "\t")
np.savetxt("5c_" + strain + "_hist_bin_100.txt", histBinMatrix, fmt = "%s", delimiter = "\t")
'''


#########################################################################################################
# Sort into 3 classes: one cluster, two clusters, multiple clusters
#########################################################################################################

'''
# Change to kp as required
strain = "ec"
histBin = np.loadtxt("s2a_classification_results/" + "5c_" + strain + "_hist_bin_100.txt", dtype = str)
silScore = np.loadtxt("s2a_classification_results/" + "5a_" + strain + "_sil_score.txt", dtype = str)

print histBin
print silScore

outputHeader = histBin[0, 0:4]
outputHeader = np.append(outputHeader, "type")
outputMatrix = histBin[1:, 0:4]

silScoreVals = silScore[1:, 4].astype(float)

numSamples = outputMatrix.shape[0]

afMatrix = histBin[1:, 4:].astype(float)
classifyVector = np.zeros(numSamples, dtype = object)

for sample in range(numSamples):

	# Assign samples with >= 90% of SNPs at or above 0.9 AF to "one strain"
	if np.sum(afMatrix[sample, 90:]) / np.sum(afMatrix[sample, :]) >= 0.90:
		classifyVector[sample] = "one_class"

	# Of those not assigned to "one strain", assign samples with silhouette score >= 0.8 for k = 2 to "two strains"
	elif silScoreVals[sample] >= 0.8:
		classifyVector[sample] = "two_class"

	# Assign all other samples to "multiple strains"	
	else:
		classifyVector[sample] = "multiple_class"

outputMatrix = np.hstack((outputMatrix, classifyVector.reshape(-1, 1)))
outputMatrix = np.vstack((outputHeader, outputMatrix))

np.savetxt("5d_" + strain + "_classify.txt", outputMatrix, fmt = "%s", delimiter = "\t")
'''

#########################################################################################################
# Count transitions
#########################################################################################################

'''
# Change to Kleb as required
strain = "Ecoli"

classifyMatrix = np.loadtxt("s2a_classification_results/updated_old_values/a2_" + strain + "_classify_grey_split.txt", skiprows = 1, dtype = str)

print (classifyMatrix.shape)

transAll = np.zeros((5, 5), dtype = int)
transSub = np.zeros((5, 5), dtype = int)
transFamMem = np.zeros((5, 5), dtype = int)

strayCount = 0

for row in classifyMatrix:
	for pos in range(2, 13):

		if row[pos] == "G" and row[pos + 1] == "G":
			entryIndex = (0, 0)
		elif (row[pos] == "G" and row[pos + 1] == "GM"):
			entryIndex = (0, 1)
		elif (row[pos] == "G" and row[pos + 1] == "B"):
			entryIndex = (0, 2)
		elif row[pos] == "G" and row[pos + 1] == "Y":
			entryIndex = (0, 3)
		elif row[pos] == "G" and row[pos + 1] == "R":
			entryIndex = (0, 4)

		elif row[pos] == "GM" and row[pos + 1] == "G":
			entryIndex = (1, 0)
		elif (row[pos] == "GM" and row[pos + 1] == "GM"):
			entryIndex = (1, 1)
		elif (row[pos] == "GM" and row[pos + 1] == "B"):
			entryIndex = (1, 2)
		elif row[pos] == "GM" and row[pos + 1] == "Y":
			entryIndex = (1, 3)
		elif row[pos] == "GM" and row[pos + 1] == "R":
			entryIndex = (1, 4)

		elif row[pos] == "B" and row[pos + 1] == "G":
			entryIndex = (2, 0)
		elif row[pos] == "B" and row[pos + 1] == "GM":
			entryIndex = (2, 1)
		elif row[pos] == "B" and row[pos + 1] == "B":
			entryIndex = (2, 2)
		elif row[pos] == "B" and row[pos + 1] == "Y":
			entryIndex = (2, 3)
		elif row[pos] == "B" and row[pos + 1] == "R":
			entryIndex = (2, 4)

		elif row[pos] == "Y" and row[pos + 1] == "G":
			entryIndex = (3, 0)
		elif row[pos] == "Y" and row[pos + 1] == "GM":
			entryIndex = (3, 1)
		elif row[pos] == "Y" and row[pos + 1] == "B":
			entryIndex = (3, 2)
		elif row[pos] == "Y" and row[pos + 1] == "Y":
			entryIndex = (3, 3)
		elif row[pos] == "Y" and row[pos + 1] == "R":
			entryIndex = (3, 4)

		elif row[pos] == "R" and row[pos + 1] == "G":
			entryIndex = (4, 0)
		elif row[pos] == "R" and row[pos + 1] == "GM":
			entryIndex = (4, 1)
		elif row[pos] == "R" and row[pos + 1] == "B":
			entryIndex = (4, 2)
		elif row[pos] == "R" and row[pos + 1] == "Y":
			entryIndex = (4, 3)
		elif row[pos] == "R" and row[pos + 1] == "R":
			entryIndex = (4, 4)

		else:
			strayCount += 1
			continue

		transAll[entryIndex[0], entryIndex[1]] += 1

		if row[1] == "T":
			transSub[entryIndex[0], entryIndex[1]] += 1
		else:
			transFamMem[entryIndex[0], entryIndex[1]] += 1

print (strayCount)
print (transAll)
print (transSub)
print (transFamMem)
'''

#########################################################################################################
# Generate string for markovchain in R
#########################################################################################################

'''
# Change to Kleb as required
strain = "Ecoli"

# Change to subj as required
sampleType = "fam"

statusMatrix = np.loadtxt("s2a_classification_results/updated_old_values/a3_" + strain + "_classify_grey_split_" + sampleType + ".txt", dtype = str)

transitions = statusMatrix[1:, 2:]
seqString = ""
for row in transitions:
	for state in row:
		if state == "N":
			seqString = seqString + "NA, "
		elif state == "G":
			seqString = seqString + "NA, "
		elif state == "GM":
			seqString = seqString + "\"missing\", "
		elif state == "B":
			seqString = seqString + "\"one\", "
		elif state == "Y":
			seqString = seqString + "\"two\", "
		elif state == "R":
			seqString = seqString + "\"multiple\", "
	seqString = seqString + "NA, "

print seqString
'''



