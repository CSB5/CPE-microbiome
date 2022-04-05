import numpy as np
import os

cutoff = 0.1

abun = np.loadtxt("Enterobact_merged.dat", dtype = str)
metaData = np.loadtxt("Supplementary_file_1.txt", usecols = (1,4,7), dtype = str)

abunFileIDs = abun[0, 1:]
metaIDs = metaData[1:, 2]
commonIDs = np.intersect1d(abunFileIDs, metaIDs)

outputTable = np.zeros((1,6), dtype = object)
outputTable[0, 0:7] = ["SubjectType", "IndivCode", "SequencingLibraryID", "EcRelAbun", "KpRelAbun", "EnteroRelAbun"]

# Extract relevant abundance values
for row in metaData:

	if row[2] not in commonIDs:
		continue

	outputRow = np.zeros((1,6), dtype = object)
	outputRow[0, 0:3] = row
	abunFileCol = np.where(abun[0, :] == row[2])[0][0]
	outputRow[0, 3] = abun[2, abunFileCol].astype(float)
	outputRow[0, 4] = abun[4, abunFileCol].astype(float)
	outputRow[0, 5] = np.sum(abun[1:, abunFileCol].astype(float))

	outputTable = np.vstack((outputTable, outputRow))

np.savetxt("relAbun.txt", outputTable, fmt = "%s", delimiter = "\t")

outputStDevTable = np.zeros((1,5), dtype = object)
outputStDevTable[0, 0:6] = ["SubjectType", "IndivCode", "EcAbunStDev", "KpAbunStDev", "EnteroAbunStDev"]

indivList = np.unique(outputTable[1:, 1])

for indiv in indivList:
	extractedRows = outputTable[outputTable[:, 1] == indiv, :]

	# Ignore subjects with only one timepoint
	if extractedRows.shape[0] == 1:
		continue

	ecVals = extractedRows[:, 3]
	ecVals = ecVals[ecVals > cutoff]

	kpVals = extractedRows[:, 4]
	kpVals = kpVals[kpVals > cutoff]

	enteroVals = extractedRows[:, 5]
	enteroVals = enteroVals[enteroVals > cutoff]

	outputStDevRow = np.zeros((1,5), dtype = object)
	outputStDevRow[0, 0:2] = extractedRows[0, 0:2]

	# Ignore subjects with only one timepoint after filtering for the minimum abundance value
	if len(ecVals) >= 2:
		outputStDevRow[0, 2] = np.std(ecVals, ddof = 1)
	else:
		outputStDevRow[0, 2] = "NA"

	if len(kpVals) >= 2:
		outputStDevRow[0, 3] = np.std(kpVals, ddof = 1)
	else:
		outputStDevRow[0, 3] = "NA"

	if len(enteroVals) >= 2:
		outputStDevRow[0, 4] = np.std(enteroVals, ddof = 1)
	else:
		outputStDevRow[0, 4] = "NA"

	outputStDevTable = np.vstack((outputStDevTable, outputStDevRow))

np.savetxt("relAbunStDev.txt", outputStDevTable, fmt = "%s", delimiter = "\t")


















