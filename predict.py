from data import *
import scipy.ndimage as spn


class predictDose(object):
    def __init__(self,dat,refDose):
        if isinstance(dat,data):
            self.data = dat
        else:
            print 'data type wrong'
        self.refDose = refDose
        self.predictedDose = np.ones(self.data.nVox)
        self.structureWeights = {}
        for s in range(self.data.nStructures):
            self.structureWeights[self.data.structureNames[s]] = (0.,0.)

        self.voxelThresh = self.data.thresh.copy()
        self.voxelOverPenalty = self.data.overPenalty.copy()
        self.voxelUnderPenalty = self.data.underPenalty.copy()
        self.sortedStructureIndices = []
        self.refDosePerStruct = []

    def genDosePerStruct(self,source,destination):
        del destination[:]
        for s in range(self.data.nStructures):
            destination.append(source[self.data.structureVoxels[s]].copy())

    def genSortedIndicesPerStruct(self,arraysToSort,destinationForIndices, sortDirection=1):
        del destinationForIndices[:]
        for s in range(self.data.nStructures):
            destinationForIndices.append(np.argsort(arraysToSort[s], axis=0)[::sortDirection])

    def genSortedDosePerStruct(self, doseToSort, destination, sortDirection = 1):
        del destination[:]
        for s in range(self.data.nStructures):
            destination.append(np.sort(doseToSort[self.data.structureVoxels[s]])[::sortDirection])

    def predictDose(self):
        print 'write dose prediction function'

    def genStructureWeightingFromArea(self,oarExponent=1., ptvOverdoseExponent=1., ptvUnderdoseExponent=1.):
        print 'generating structure weights based on ref dose DVH area'
        for s in range(self.data.nStructures):
            wO, wU = 0., 0.
            if self.data.structureNames[s] not in set(self.data.PTVNames) and self.data.structureVoxels[s].size >0:
                wO += 1.+np.average(self.refDose[self.data.structureVoxels[s]]) ** oarExponent
            elif self.data.structureVoxels[s].size>0:
                wO += 1.+np.average(self.refDose[self.data.structureVoxels[s]]) ** ptvOverdoseExponent
                wU += 1.+np.average(self.refDose[self.data.structureVoxels[s]]) ** ptvUnderdoseExponent
            self.structureWeights[self.data.structureNames[s]] = (wU, wO)

    def setStructureWeightingFromValues(self, OARU=-1, OARO=-1, PTVU=-1, PTVO=-1,oarExponent=1., ptvOverdoseExponent=1., ptvUnderdoseExponent=1.):
        for s in range(self.data.nStructures):
            wO, wU = 0., 0.
            if self.data.structureNames[s] not in set(self.data.PTVNames) and self.data.structureVoxels[s].size > 0:
                wO = OARO ** oarExponent
            elif self.data.structureVoxels[s].size>0:
                wO = PTVO ** ptvOverdoseExponent
                wU = PTVU ** ptvUnderdoseExponent
            self.structureWeights[self.data.structureNames[s]] = (wU, wO)

    def setStructureWeightingFromDict(self,specificWeights={}):
        for s in range(self.data.nStructures):
            if self.data.structureNames[s] in specificWeights:
                self.structureWeights[s] = specificWeights[self.data.structureNames[s]]

    def updateThreshAndWeights(self):
        print 'write thresh update function'

    def updateThreshAndWeights(self, sortedRefDose, voxelRedistributionIndices, updateTargetDose=False, targetScalingFactor=0.5):
        print 'Updating thresh and weights'
        for s in range(self.data.nStructures):
            sReal = self.data.structureNamesInverse[self.data.weightPriorityDict[s]]

            if (self.data.structureNames[sReal] not in set(self.data.PTVNames)) and self.data.structureVoxels[sReal].size > 0:
                print 'OAR: ',self.data.structureNames[sReal]
                self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = sortedRefDose[sReal]

                self.voxelUnderPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / \
                                                                                                    self.data.structureVoxels[
                                                                                                            sReal].size
                self.voxelOverPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                    sReal].size

            elif self.data.structureVoxels[sReal].size>0:
                #keeps PTV the same!
                if updateTargetDose:
                    #self.voxelThresh[self.data.structureVoxels[sReal]] = np.minimum(sortedRefDose[sReal],self.data.threshDict[self.data.weightPriorityDict[s]]) + targetScalingFactor * np.abs(self.data.threshDict[self.data.weightPriorityDict[s]] - sortedRefDose[sReal])
                    self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = (1.-targetScalingFactor) * sortedRefDose[sReal] + 1. * targetScalingFactor * self.data.threshDict[self.data.weightPriorityDict[s]]
                    self.voxelUnderPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / self.data.structureVoxels[
                        sReal].size
                    self.voxelOverPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                        sReal].size
                else:
                    self.voxelThresh[self.data.structureVoxels[sReal]] = self.data.threshDict[self.data.weightPriorityDict[s]]
                    self.voxelUnderPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / self.data.structureVoxels[
                        sReal].size
                    self.voxelOverPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                        sReal].size

        self.data.thresh = self.voxelThresh.copy()
        self.data.underPenalty = self.voxelUnderPenalty.copy()
        self.data.overPenalty = self.voxelOverPenalty.copy()

class predictDoseDistance(predictDose):
    def __init__(self,dat,refDose):
        super(predictDoseDistance,self).__init__(dat,refDose)

        self.xDim = self.data.voxDim[0]
        self.yDim = self.data.voxDim[1]
        self.zDim = self.data.voxDim[2]
        self.distancesGrid = None
        self.distancesFull = np.zeros(self.data.nVoxFull)
        self.distances = np.zeros(self.data.nVox)
        self.distanceMask = np.ones(self.data.nVoxFull)
        self.structureDistancesSortedIndices = []
        self.structureDistances = []
        self.structureDoseSorted = []

    def predictDose(self):
        self.buildMask()
        self.calcAndSortDistances()

        self.genDosePerStruct(self.distances, self.structureDistances)
        self.genSortedIndicesPerStruct(self.structureDistances, self.structureDistancesSortedIndices,sortDirection = 1)
        self.genSortedDosePerStruct(self.refDose, self.structureDoseSorted,sortDirection=-1)
        self.updateThreshAndWeights(self.structureDoseSorted,self.structureDistancesSortedIndices)

    def buildMask(self):
        self.distanceMask = np.ones(self.data.nVoxFull)
        for pname in self.data.PTVNames:
            self.distanceMask[self.data.structureVoxelsFull[self.data.structureNamesInverse[pname]]] = 0
        self.distanceMask = self.distanceMask.reshape(self.data.voxDim, order='F')

    def calcAndSortDistances(self):
        self.distancesGrid = spn.distance_transform_edt(self.distanceMask, sampling=self.data.voxSamplingSizes[:])
        self.distancesFull = self.distancesGrid.flatten(order='F').copy()
        self.distances = self.distancesFull[self.data.voxelIndicesInFull]

    # def updateThreshAndWeights(self,updateTargetDose=False, targetScalingFactor=0.5):
    #     print 'Updating thresh and weights'
    #     for s in range(self.data.nStructures):
    #         sReal = self.data.structureNamesInverse[self.data.weightPriorityDict[s]]
    #
    #         if (self.data.structureNames[sReal] not in set(self.data.PTVNames)) and self.data.structureVoxels[sReal].size > 0:
    #             print 'OAR: ',self.data.structureNames[sReal]
    #             self.voxelThresh[self.data.structureVoxels[sReal][self.structureDistancesSortedIndices[sReal]]] = self.structureDoseSorted[sReal]
    #
    #             self.voxelUnderPenalty[self.data.structureVoxels[sReal][self.structureDistancesSortedIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / \
    #                                                                                                 self.data.structureVoxels[
    #                                                                                                         sReal].size
    #             self.voxelOverPenalty[self.data.structureVoxels[sReal][self.structureDistancesSortedIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
    #                 sReal].size
    #
    #         elif self.data.structureVoxels[sReal].size>0:
    #             #keeps PTV the same!
    #             self.voxelThresh[self.data.structureVoxels[sReal]] = self.data.threshDict[self.data.weightPriorityDict[s]]
    #             self.voxelUnderPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / self.data.structureVoxels[
    #                 sReal].size
    #             self.voxelOverPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
    #                 sReal].size
    #
    #     self.data.thresh = self.voxelThresh.copy()
    #     self.data.underPenalty = self.voxelUnderPenalty.copy()
    #     self.data.overPenalty = self.voxelOverPenalty.copy()