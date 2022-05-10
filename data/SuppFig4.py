import numpy as np
import os

abun = np.loadtxt("SuppFig4_NoEntero_Profiles.tsv", dtype = str)

subjects = abun[0, 1:]
genera = abun[1:, 0]
abunFreq = abun[1:, 1:].astype(float)

shannonScoreVec = np.zeros(len(subjects), dtype = float)
countsVec = np.zeros(len(subjects), dtype = int)

for i in range(len(subjects)):
	abunVec = abunFreq[:, i]

	abunVecProp = abunVec / 100

	# Normalize to 1
	nonEnteroProp = np.sum(abunVecProp)
	abunVecProp = abunVecProp / nonEnteroProp

	# Calculate Shannon diversity
	with np.errstate(divide = 'ignore'):
		abunVecLn = np.log(abunVecProp)
	entropyVec = np.nan_to_num(abunVecLn) * abunVecProp
	shannonScoreVec[i] = np.sum(-entropyVec)

	# Calculate richness
	abunVecCount = abunVec > 0
	countsVec[i] = np.sum(abunVecCount)

outputMatrix = np.vstack((subjects, shannonScoreVec, countsVec)).T
np.savetxt("SuppFig4_NoEntero_Summary", outputMatrix, delimiter = "\t", fmt = "%s")