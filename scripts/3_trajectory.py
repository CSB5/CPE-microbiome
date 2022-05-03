import numpy as np
import os
from distutils.dir_util import copy_tree
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

def listdir_nohidden(path):
    return [f for f in sorted(os.listdir(path)) if not f.startswith('.')]

############################################
# Produce trajectory of the union of SNPs
############################################

'''
# Change to Kleb as required
strain = "Ecoli"

indivList = listdir_nohidden(strain)

for indiv in indivList:

	# Get all timepoints across a single individual
	tpList = listdir_nohidden(strain + "/" + indiv + "/metagenomes_jon/")

	print (tpList)

	count = 0
	for tp in tpList:

		snpMatrix = np.loadtxt(strain + "/" + indiv + "/metagenomes_jon/" + tp, usecols = (2,7), dtype = object)

		vals, inverse, countUnique = np.unique(snpMatrix[:, 0].astype(int), return_inverse = True, return_counts = True)
		idx_vals_unique = np.where(countUnique == 1)[0]
		vals_unique = vals[idx_vals_unique]
		snpMatrix = snpMatrix[np.isin(snpMatrix[:, 0].astype(int), vals_unique), :]

		if count == 0:
			uniqueSnpList = snpMatrix[:, 0].astype(int)
		else:
			uniqueSnpList = np.unique(np.concatenate((uniqueSnpList, snpMatrix[:, 0].astype(int))))

		count += 1

	trajectMatrix = np.zeros((len(uniqueSnpList) + 1, len(tpList) + 1), dtype = object)
	trajectMatrix[0, 0] = "traject"
	trajectMatrix[1:, 0] = uniqueSnpList
	trajectMatrix[0, 1:] = [tp.split("_")[0] for tp in tpList]

	print (trajectMatrix)

	count = 0
	for tp in tpList:

		snpMatrix = np.loadtxt(strain + "/" + indiv + "/metagenomes_jon/" + tp, usecols = (2,7), dtype = object)

		vals, inverse, countUnique = np.unique(snpMatrix[:, 0].astype(int), return_inverse = True, return_counts = True)
		idx_vals_unique = np.where(countUnique == 1)[0]
		vals_unique = vals[idx_vals_unique]
		snpMatrix = snpMatrix[np.isin(snpMatrix[:, 0].astype(int), vals_unique), :]

		afVector = [row[1] for row in snpMatrix]
		indexVec = np.isin(uniqueSnpList, snpMatrix[:, 0].astype(int))

		trajectMatrix[1:, count + 1][indexVec] = afVector

		count +=1

	print (trajectMatrix)

	np.savetxt(strain + "/" + indiv + "/trajectUnion.txt", trajectMatrix, fmt = "%s", delimiter = "\t")
'''

############################################################################
# Trajectory analysis (intersection of SNPs in an individual + 0.9 filter)
############################################################################

'''
# Change to Kleb as required
strain = "Ecoli"

indivList = listdir_nohidden(strain)

for indiv in indivList:

	trajectFile = np.loadtxt(strain + "/" + indiv + "/trajectUnion.txt", dtype = str)

	print (trajectFile)

	# Boolean to enforce presence at all timepoints
	boolVect = np.zeros(trajectFile.shape[0], dtype = bool)
	boolVect[0] = True

	count = 1
	for row in trajectFile[1:, :]:
		if "0" not in row:
			boolVect[count] = True
		count += 1

	newTrajectFile = trajectFile[boolVect, :]

	varyBool = np.zeros(newTrajectFile.shape[0], dtype = bool)
	varyBool[0] = True
	# Exclude SNPs with mean AF >= 0.9 across all available timepoints
	varyBool[1:] = np.mean(newTrajectFile[1:, 1:].astype(float), axis = 1) < 0.9

	newTrajectFile = newTrajectFile[varyBool, :]

	np.savetxt(strain + "/" + indiv + "/trajectIntsctVary.txt", newTrajectFile, fmt = "%s", delimiter = "\t")
'''

#########################################################################################################
# Add labels: Which cluster, and whether present in isolate
#########################################################################################################

'''
# Change to Kleb as required
strain = "Ecoli"

# Change as required
indiv = "1674t"

numSamp = path = strain + "/" + indiv + "/"

isoFileList = listdir_nohidden(path + "isolates/")

# Get union of all SNPs across the isolates at different timepoints
count = 0
for isoFile in isoFileList:
	isoFileSnp = np.loadtxt(path + "isolates/" + isoFile, usecols = (1,), skiprows = 1, dtype = int)
	if count == 0:
		isolateSnpUnion = isoFileSnp
	else:
		isolateSnpUnion = np.unique(np.concatenate((isolateSnpUnion, isoFileSnp)))
	count += 1

isolateSnpUnion = np.sort(isolateSnpUnion)

trajectFile = np.loadtxt(path + "trajectIntsctVary.txt", dtype = object)
dim = trajectFile.shape[1] - 1

# Parameters for DBSCAN
epsVal = 0.2
minSampVal = 2 * dim

rowHead = trajectFile[0, :]
rowHead = np.append(rowHead, ["class", "in_isolate"])
colHead = np.reshape(trajectFile[1:, 0], (-1, 1)).astype(int)
trajectFile = trajectFile[1:, 1:].astype(float)

inIso = np.isin(colHead, isolateSnpUnion) * 1

dbClust = DBSCAN(eps = epsVal, min_samples = minSampVal).fit(trajectFile)
# Noisy samples labeled -1, add 2 across the board to start index from 1
labels = np.reshape(dbClust.labels_ + 2, (-1, 1)).astype(object)

print (labels.shape)

outputFile = np.hstack((colHead, trajectFile.astype(object), labels, inIso))
outputFile = np.vstack((rowHead, outputFile))

print outputFile

np.savetxt(path + "trajectIntsctVary_labels.txt", outputFile, fmt = "%s", delimiter = "\t")
'''
	