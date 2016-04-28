import os
import numpy as np
import scipy.sparse as sps
import scipy.io as sio


class data(object):
    def __init__(self, dataDict):
        # read in shared data bits
        print 'building common data'
        self.workingDir = dataDict['workingDir']
        self.dataDir = dataDict['dataDir']
        self.beamNumbers = dataDict['beams'][:]
        self.nBeams = dataDict['nbeams']
        self.voxSamplingSizes = dataDict['voxSampling'][:]
        self.voxDim = dataDict['voxDim']
        self.PTVNames = dataDict['ptvnames']
        self.aperLimit = dataDict['aperLimitFlag'] or self.nBeams
        self.useGPU = dataDict['useGPU']
        self.modelType = dataDict['modelType']
        self.dualThresh = dataDict['dualThresh']

        self.structureNames = []
        self.structureVoxelsFull = []
        self.nBeamlets = 0
        self.nBPB = []

        voiFilenameList = []
        for file in os.listdir(self.dataDir):
            if file.endswith('VOILIST.mat'):
                voiFilenameList.append(file)
        self.nStructures = len(voiFilenameList)

        initialMat = sio.loadmat(self.dataDir + 'Gantry' + str(self.beamNumbers[0]) + '_Couch0_D.mat')
        self.nVoxFull = initialMat['D'].shape[0]

        onesWhereVoxels = np.zeros(self.nVoxFull, dtype='int64')
        for s in range(len(voiFilenameList)):
            sMat = sio.loadmat(self.dataDir + voiFilenameList[s])
            self.structureNames.append(voiFilenameList[s][:-12])
            self.structureVoxelsFull.append(sMat['v'] - 1)
            onesWhereVoxels[self.structureVoxelsFull[s]] = 1

        self.nVox = int(onesWhereVoxels.sum())
        self.voxelIndicesInFull = onesWhereVoxels.nonzero()[0]
        self.voxelIndicesInTrunc = np.zeros(self.nVoxFull, dtype='int64')
        for i in range(self.nVox):
            self.voxelIndicesInTrunc[self.voxelIndicesInFull[i]] = i

        # read in beams
        self.dijList = []
        for b in self.beamNumbers:
            print 'reading in beam', b
            dMat = sio.loadmat(self.dataDir + 'Gantry' + str(b) + '_Couch0_D.mat')
            d = sps.csc_matrix(dMat['D'][self.voxelIndicesInFull, :])
            self.nBeamlets += d.shape[1]
            self.nBPB.append(d.shape[1])
            self.dijList.append(d)

        self.nVox = self.dijList[0].shape[0]
        self.nBPBcum = [0] + [i for i in np.cumsum(self.nBPB)]

        print 'building objective'
        self.weightUODict = dataDict['weightUnderOverDict']
        self.weightPriorityDict = dataDict['weightPriorityDict']
        self.weightPriorityDictInverse = {}
        for key, value in self.weightPriorityDict.iteritems():
            self.weightPriorityDictInverse[value] = key
        self.threshDict = dataDict['threshDict']
        self.structureNamesInverse = {}
        for s in range(self.nStructures):
            self.structureNamesInverse[self.structureNames[s]] = s

        self.structureVoxelsOverlap = []
        for s in range(self.nStructures):
            self.structureVoxelsOverlap.append(self.voxelIndicesInTrunc[self.structureVoxelsFull[s]])

        self.singleStructureAssignmentVector = np.zeros(self.nVox) - 1
        for s in range(self.nStructures):
            sReal = self.structureNamesInverse[self.weightPriorityDict[s]]
            self.singleStructureAssignmentVector[self.structureVoxelsOverlap[sReal]] = sReal

        self.structureVoxels = []
        for s in range(self.nStructures):
            self.structureVoxels.append(np.where(self.singleStructureAssignmentVector == s)[0])

        for s in range(self.nStructures):
            print s, self.structureVoxels[s].size, self.structureVoxelsOverlap[s].size, self.structureNames[s]

        self.underPenalty = np.zeros(self.nVox, dtype='float64')
        self.overPenalty = np.zeros(self.nVox, dtype='float64')
        self.thresh = np.zeros(self.nVox, dtype='float64')

        for s in range(self.nStructures):
            if self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]].size > 0:
                self.underPenalty[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = 1. * \
                                                                                                                  self.weightUODict[
                                                                                                                      self.weightPriorityDict[
                                                                                                                          s]][
                                                                                                                      0] / \
                                                                                                                  self.structureVoxels[
                                                                                                                      self.structureNamesInverse[
                                                                                                                          self.weightPriorityDict[
                                                                                                                              s]]].size
                self.overPenalty[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = 1. * \
                                                                                                                 self.weightUODict[
                                                                                                                     self.weightPriorityDict[
                                                                                                                         s]][
                                                                                                                     1] / \
                                                                                                                 self.structureVoxels[
                                                                                                                     self.structureNamesInverse[
                                                                                                                         self.weightPriorityDict[
                                                                                                                             s]]].size
                self.thresh[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = \
                    self.threshDict[self.weightPriorityDict[s]]

        if dataDict['threshOverride'] is not None:
            self.thresh = dataDict['threshOverride']
        if dataDict['overWOverride'] is not None:
            self.overPenalty = dataDict['overWOverride']
        if dataDict['underWOverride'] is not None:
            self.underPenalty = dataDict['underWOverride']
        if dataDict['comparisonDose'] is not None:
            self.comparisonDose = dataDict['comparisonDose']
        else:
            self.comparisonDose = None

        if 'underTOverride' not in dataDict.keys():
            self.underThresh = self.thresh.copy()
        else:
            self.underThresh = dataDict['underTOverride']
        if 'overTOverride' not in dataDict.keys():
            self.overThresh = self.thresh.copy()
        else:
            self.overThresh = dataDict['overTOverride']

        # set base penalty
        if 'basePenalty' not in dataDict.keys():
            self.basePenalty = False
        else:
            self.basePenalty = dataDict['basePenalty']

        if 'basePenaltyWeightOAR' not in dataDict.keys():
            self.basePenaltyWeightOAR = 1.
        else:
            self.basePenaltyWeightOAR = dataDict['basePenaltyWeightOAR']

        if 'basePenaltyWeightPTV' not in dataDict.keys():
            self.basePenaltyWeightPTV = 1.
        else:
            self.basePenaltyWeightPTV = dataDict['basePenaltyWeightPTV']

        # build base penalty threshold
        self.buildBasePenaltyVectors()

        # read in modality part
        self.readInModalityData(dataDict)

    def buildBasePenaltyVectors(self, oarWeight=None, ptvWeight=None):
        if oarWeight is None:
            oarWeight = self.basePenaltyWeightOAR
        if ptvWeight is None:
            ptvWeight = self.basePenaltyWeightPTV
        self.baseThresh = np.zeros(self.nVox, dtype='float64')
        self.basePenaltyOver = np.zeros(self.nVox, dtype='float64')
        self.basePenaltyUnder = np.zeros(self.nVox, dtype='float64')
        for s in range(self.nStructures):
            if self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]].size > 0:
                # IF STRUCTURE IS OAR THEN USE OAR PENALTY, O/W USE TARGET
                over, under = oarWeight, oarWeight
                if self.structureNamesInverse[self.weightPriorityDict[s]] in self.PTVNames:
                    over, under = ptvWeight, ptvWeight

                self.basePenaltyUnder[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = 1. * \
                                                                                                                  under / \
                                                                                                                  self.structureVoxels[
                                                                                                                      self.structureNamesInverse[
                                                                                                                          self.weightPriorityDict[
                                                                                                                              s]]].size
                self.basePenaltyOver[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = 1. * \
                                                                                                                  over / \
                                                                                                                  self.structureVoxels[
                                                                                                                      self.structureNamesInverse[
                                                                                                                          self.weightPriorityDict[
                                                                                                                              s]]].size

                self.baseThresh[self.structureVoxels[self.structureNamesInverse[self.weightPriorityDict[s]]]] = self.threshDict[self.weightPriorityDict[s]]


    def updateOverUnderThresh(self, newThreshU=None, newThreshO=None):
        if newThreshU is None:
            newThreshU = self.thresh.copy()
        if newThreshO is None:
            newThreshO = self.thresh.copy()
        self.overThresh = newThreshO.copy()
        self.underThresh = newThreshU.copy()

    def readInModalityData(self, dataDict):
        print 'please write child function'

    def updateObjParamsFromFile(self, filename=None):
        if filename is None:
            print 'you did this wrong'
        else:
            newData = sio.loadmat(filename)
            self.thresh = newData['thresh'].flatten()
            self.underPenalty = newData['under'].flatten()
            self.overPenalty = newData['over'].flatten()
            print 'Objective parameters updated from', filename


class data_fmo(data):
    def __init__(self, dataDict):
        if dataDict['modelType'] != 'fmo':
            print 'ERROR IN MODEL TYPE'
            exit()
        super(data_fmo, self).__init__(dataDict)

    def readInModalityData(self, dataDict):
        print 'child Modality data readin'


class data_dao(data):
    def __init__(self, dataDict):

        if dataDict['modelType'] != 'vmat' and dataDict['modelType'] != 'dao':
            print 'ERROR IN MODEL TYPE'
            exit()
        super(data_dao, self).__init__(dataDict)

        print 'calling vmat-specific data constructor'
        self.rowStep = dataDict['vmatRowStep']
        self.colStep = dataDict['vmatColStep']
        self.leafDistPerDeg = dataDict['leafDistPerDeg']
        self.degPerStep = 360. / self.nBeams
        self.vmatFlag = dataDict['vmatFlag']

        # extract row/column data and beam data
        self.rowStartPB = []
        self.rowEndPB = []
        self.numRows = []

        print 'getting row col max and min dimensions'
        minRow = 1000
        maxRow = -1
        minCol = 1000
        maxCol = -1
        for b in self.beamNumbers:
            print 'Scanning in beam', b
            # beamlet positions
            rowMat = sio.loadmat(self.dataDir + 'Gantry' + str(b) + '_Couch0_BEAMINFO.mat')

            beamx = rowMat['x'].flatten()
            beamy = rowMat['y'].flatten()
            minRow = min(minRow, np.amin(beamx))
            maxRow = max(maxRow, np.amax(beamx))
            minCol = min(minCol, np.amin(beamx))
            maxCol = max(maxCol, np.amax(beamx))
        spacedColMin = 0
        spacedColMax = int((maxCol - minCol) * self.colStep) + 1
        spacedRowMin = 0
        spacedRowMax = int((maxRow - minRow) * self.rowStep) + 1
        print 'beamlet spacing:', spacedColMin, spacedColMax, spacedRowMin, spacedRowMax
        self.spacedColMax = spacedColMax
        self.spacedColMin = spacedColMin
        self.spacedRowMin = spacedRowMin
        self.spacedRowMax = spacedRowMax
        self.DkjRMP = None
        self.Dkj = None

        print 'reading in beamlet positions'
        self.rowIndToActual = []
        self.rowActualToInd = []
        self.colStartPerRow = []
        for b in self.beamNumbers:
            # beamlet positions
            rowMat = sio.loadmat(self.dataDir + 'Gantry' + str(b) + '_Couch0_BEAMINFO.mat')
            beamx = rowMat['x'].flatten()
            beamy = rowMat['y'].flatten()
            beamxunique = np.unique(beamx)
            numrows = beamxunique.size
            rs = np.zeros(numrows, dtype='int32')
            re = np.zeros(numrows, dtype='int32')

            rIndToActual = np.zeros(numrows, dtype='int32')
            rActualToInd = np.zeros(spacedRowMax, dtype='int32') - 1
            rColStartPerRow = np.zeros(numrows, dtype='int32')

            # todo extract beam position in this loop
            for u in range(numrows):
                indices = np.where(beamx == beamxunique[u])[0]
                rs[u] = indices.min()
                re[u] = indices.max() + 1
                rIndVal = int((beamxunique[u] - minRow) * self.rowStep)
                rActualToInd[rIndVal] = u
                rIndToActual[u] = rIndVal
                rColStartPerRow[u] = int((beamy[indices.min()] - minCol) * self.colStep)
            self.rowIndToActual.append(rIndToActual)
            self.rowActualToInd.append(rActualToInd)
            self.colStartPerRow.append(rColStartPerRow)
            self.rowStartPB.append(rs)
            self.rowEndPB.append(re)
            self.numRows.append(numrows)

    def readInModalityData(self, dataDict):
        print 'child Modality data readin VMAT'
