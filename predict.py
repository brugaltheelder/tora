from data import *
import operator
import scipy.ndimage as spn
from npgamma import *
import matplotlib.pyplot as plt








def genAlphas(dat,PTVBase = 1000., OARBase = 10.,modifier=['default']):
    # assign bases
    assert(isinstance(dat,data_fmo))



    for n in range(dat.nStructures):
        if dat.structureNames[n] in dat.PTVNames:
            dat.weightUODict[dat.structureNames[n]] = (4*PTVBase,PTVBase)
        else:
            dat.weightUODict[dat.structureNames[n]] = (1.,OARBase)

    # assign modifiers

    if 'PTV' in modifier:
        for s in dat.PTVNames:
            dat.weightUODict[s] = (dat.weightUODict[s][0]**2,dat.weightUODict[s][1]**2)
    if 'Bladder' in modifier:
        dat.weightUODict['Bladder'] = (dat.weightUODict['Bladder'][0]**2,dat.weightUODict['Bladder'][1]**2)
    if 'Rectum' in modifier:
        dat.weightUODict['Rectum'] = (dat.weightUODict['Rectum'][0]**2,dat.weightUODict['Rectum'][1]**2)

    dat.buildPenaltyVectors(dat.underPenalty,dat.overPenalty,dat.thresh)




class predictDose(object):
    def __init__(self,dat,refDose):
        if isinstance(dat, data):
            self.data = dat
        else:
            print 'data type wrong'
        self.refDose = refDose.copy()
        self.predictedDose = np.ones(self.data.nVox)
        self.structureWeights = {}
        for s in range(self.data.nStructures):
            self.structureWeights[self.data.structureNames[s]] = (0.,0.)

        self.voxelThresh = self.data.thresh.copy()
        self.voxelOverPenalty = self.data.overPenalty.copy()
        self.voxelUnderPenalty = self.data.underPenalty.copy()
        self.xDim = self.data.voxDim[0]
        self.yDim = self.data.voxDim[1]
        self.zDim = self.data.voxDim[2]
        self.relPositionFull = np.zeros((self.data.nVoxFull,),dtype=[('xdim','f8'),('ydim','f8'),('zdim','f8')])
        self.relPositionFull3D = self.relPositionFull.reshape(self.data.voxDim,order='F').copy()
        for index,value in np.ndenumerate(self.relPositionFull3D):
            self.relPositionFull3D[index] = (self.xDim * index[0], self.yDim * index[1], self.zDim * index[2])
        self.relPositionFull = self.relPositionFull3D.flatten(order = 'F').copy()
        self.relPosition = self.relPositionFull[self.data.voxelIndicesInFull]


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
        #print type(doseToSort),len(doseToSort), type(destination)
        for s in range(self.data.nStructures):
            destination.append(np.sort(doseToSort[self.data.structureVoxels[s]])[::sortDirection])

    def predictDose(self):
        print 'write dose prediction function'

    def genStructureWeightingFromArea(self,sourceDose=None, oarExponent=1., ptvOverdoseExponent=1., ptvUnderdoseExponent=1., scaleDict = {}):
        print 'generating structure weights based on some dose DVH area'
        if sourceDose is not None:
            dose = sourceDose.copy()
        else:
            dose = self.refDose.copy()
        for s in range(self.data.nStructures):
            valueScalar = 1.
            if self.data.structureNames[s] in scaleDict.keys():
                valueScalar = scaleDict[self.data.structureNames[s]]
            wO, wU = 0., 0.
            if self.data.structureNames[s] not in set(self.data.PTVNames) and self.data.structureVoxels[s].size >0:
                wO += ( 1.+np.average(dose[self.data.structureVoxels[s]]) ** oarExponent ) * valueScalar
            elif self.data.structureVoxels[s].size>0:
                wO += (1.+np.average(dose[self.data.structureVoxels[s]]) ** ptvOverdoseExponent) * valueScalar
                wU += (1.+np.average(dose[self.data.structureVoxels[s]]) ** ptvUnderdoseExponent) * valueScalar
            self.structureWeights[self.data.structureNames[s]] = (wU, wO)

    def setStructureWeightingFromValues(self, OARU=-1, OARO=-1, PTVU=-1, PTVO=-1,oarExponent=1., ptvOverdoseExponent=1., ptvUnderdoseExponent=1., overrideWeightScalar=1., overrideWeightStructureName=None):
        for s in range(self.data.nStructures):
            wO, wU = 0., 0.
            if self.data.structureNames[s] not in set(self.data.PTVNames) and self.data.structureVoxels[s].size > 0:
                wO = OARO ** oarExponent
            elif self.data.structureVoxels[s].size>0:
                wO = PTVO ** ptvOverdoseExponent
                wU = PTVU ** ptvUnderdoseExponent

            if overrideWeightStructureName == self.data.structureNames[s]:
                wO = wO * overrideWeightScalar
                wU = wU * overrideWeightScalar

            self.structureWeights[self.data.structureNames[s]] = (wU, wO)

    def setStructureWeightingFromDict(self,specificWeights={}):
        for s in range(self.data.nStructures):
            if self.data.structureNames[s] in specificWeights:
                self.structureWeights[s] = specificWeights[self.data.structureNames[s]]

    def getSpatialGrid(self):
        xcoords = [1.0*i*self.data.voxSamplingSizes[0] for i in range(self.data.voxDim[0])]
        ycoords = [1.0*i*self.data.voxSamplingSizes[1] for i in range(self.data.voxDim[1])]
        zcoords = [1.0*i*self.data.voxSamplingSizes[2] for i in range(self.data.voxDim[2])]
        return np.array(xcoords), np.array(ycoords), np.array(zcoords)

    def updateThreshAndWeights(self):
        print 'write thresh update function'

    def updateThreshAndWeights(self, sortedRefDose, voxelRedistributionIndices, updateTargetDose=False, targetScalingFactor=0.5, dualThresh=True, scaleDict = {}, updateWeights=False):
        print 'Updating thresh and weights'
        voxelThreshOver = None
        voxelThreshUnder = None
        if dualThresh:
            voxelThreshOver = self.voxelThresh.copy()
            voxelThreshUnder = self.voxelThresh.copy()
        for s in range(self.data.nStructures):

            sReal = self.data.structureNamesInverse[self.data.weightPriorityDict[s]]
            valueScalar = 1.
            if self.data.weightPriorityDict[s] in scaleDict.keys():
                valueScalar = scaleDict[self.data.weightPriorityDict[s]]


            if (self.data.structureNames[sReal] not in set(self.data.PTVNames)) and self.data.structureVoxels[sReal].size > 0:
                print 'OAR: ',self.data.structureNames[sReal]
                self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = sortedRefDose[sReal]

                self.voxelUnderPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / \
                                                                                                    self.data.structureVoxels[
                                                                                                            sReal].size * valueScalar
                self.voxelOverPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                    sReal].size * valueScalar

            elif self.data.structureVoxels[sReal].size>0:
                #keeps PTV the same!
                if updateTargetDose:
                    #self.voxelThresh[self.data.structureVoxels[sReal]] = np.minimum(sortedRefDose[sReal],self.data.threshDict[self.data.weightPriorityDict[s]]) + targetScalingFactor * np.abs(self.data.threshDict[self.data.weightPriorityDict[s]] - sortedRefDose[sReal])

                    self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = (1.-targetScalingFactor) * sortedRefDose[sReal] + 1. * targetScalingFactor * self.data.threshDict[self.data.weightPriorityDict[s]]
                    if dualThresh:
                        voxelThreshOver[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = np.maximum(self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]],  self.data.threshDict[self.data.weightPriorityDict[s]])
                        voxelThreshUnder[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = np.minimum(self.voxelThresh[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]],  self.data.threshDict[self.data.weightPriorityDict[s]])
                    self.voxelUnderPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / self.data.structureVoxels[
                        sReal].size * valueScalar
                    self.voxelOverPenalty[self.data.structureVoxels[sReal][voxelRedistributionIndices[sReal]]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                        sReal].size * valueScalar
                else:
                    self.voxelThresh[self.data.structureVoxels[sReal]] = self.data.threshDict[self.data.weightPriorityDict[s]]
                    self.voxelUnderPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][0] / self.data.structureVoxels[
                        sReal].size * valueScalar
                    self.voxelOverPenalty[self.data.structureVoxels[sReal]] = 1. * self.structureWeights[self.data.structureNames[sReal]][1] / self.data.structureVoxels[
                        sReal].size * valueScalar

        self.data.thresh = self.voxelThresh.copy()
        if updateWeights:
            self.data.underPenalty = self.voxelUnderPenalty.copy()
            self.data.overPenalty = self.voxelOverPenalty.copy()
        self.data.updateOverUnderThresh(newThreshO=voxelThreshOver, newThreshU=voxelThreshUnder)



    def evalDoses(self,refDose,predDose,type='gamma',
                  gammaPrefsDict = {'distThresh':3.,'doseThreshScalar':0.03,'lowerDoseThreshScalar':0.2,'maxTestDistScalar':2.}):
        comparisonValue = -1.
        comparisonVector = np.zeros(self.data.nVox)
        if type=='gamma':
            # get metadata
            ref3D = np.zeros(self.data.nVoxFull)
            ref3D[self.data.voxelIndicesInFull] = refDose
            ref3D = ref3D.reshape(self.data.voxDim, order='F').copy()
            pred3D = np.zeros(self.data.nVoxFull)
            pred3D[self.data.voxelIndicesInFull] = predDose
            pred3D = pred3D.reshape(self.data.voxDim, order='F').copy()
            coords = self.getSpatialGrid()
            gamma = calc_gamma(coords,ref3D,coords,pred3D,gammaPrefsDict['distThresh'],
                               gammaPrefsDict['doseThreshScalar'] * np.max(ref3D),
                               lower_dose_cutoff=gammaPrefsDict['lowerDoseThreshScalar'] * np.max(pred3D),
                               distance_step_size=gammaPrefsDict['distThresh']/10.,
                               maximum_test_distance=2. * gammaPrefsDict['distThresh'])
            valid_gamma = np.zeros(self.data.nVoxFull).reshape(self.data.voxDim)
            valid_gamma[~np.isinf(gamma)] = gamma[~np.isinf(gamma)]
            return sum(valid_gamma),valid_gamma


        elif type =='squaredError':
            pass
        elif type =='sortingDistance':
            pass

        return comparisonValue, comparisonVector


    def plotSlice(self,gamma,coords,doseRef, dosePred, slice):
        max_ref_dose = np.max(doseRef)
        xc, yc, zc = coords
        cut_off_gamma = gamma.copy()
        greater_than_2_ref = (cut_off_gamma > 2.) & ~np.isinf(gamma)
        cut_off_gamma[greater_than_2_ref] = 2.
        z_i = zc[slice]



        print("Slice = {0}".format(z_i))

        plt.contourf(
            xc, yc, dosePred[:, :, slice], 30,
            vmin=0, vmax=max_ref_dose, cmap=plt.get_cmap('gist_heat'))
        plt.title("Predicted")
        plt.colorbar()
        plt.show()

        plt.contourf(
            xc, yc, doseRef[:, :, slice], 30,
            vmin=0, vmax=max_ref_dose, cmap=plt.get_cmap('gist_heat'))
        plt.title("Reference")
        plt.colorbar()
        plt.show()

        plt.contourf(
            xc, yc, cut_off_gamma[:, :, slice], 30,
            vmin=0, vmax=2, cmap=plt.get_cmap('bwr'))
        plt.title("Gamma")
        plt.colorbar()
        plt.show()








class predictDoseDistance(predictDose):
    def __init__(self,dat,refDose):
        super(predictDoseDistance, self).__init__(dat,refDose)


        self.voxelSampling = [self.xDim,self.yDim,self.zDim * self.data.zSamplingScale]
        self.distancesGrid = None
        self.distancesFull = np.zeros(self.data.nVoxFull)
        self.distances = np.zeros(self.data.nVox)
        self.distanceMask = np.ones(self.data.nVoxFull)
        self.structureDistancesSortedIndices = []
        self.structureDistances = []
        self.structureDoseSorted = []
        self.predDoseVector = None

    def predictDose(self, savePredDose=False, updateWeights = True):
        self.buildMask()
        self.calcAndSortDistances()

        self.genDosePerStruct(self.distances, self.structureDistances)
        self.genSortedIndicesPerStruct(self.structureDistances, self.structureDistancesSortedIndices,sortDirection = 1)
        self.genSortedDosePerStruct(self.refDose, self.structureDoseSorted, sortDirection=-1)
        if updateWeights:
            self.updateThreshAndWeights(self.structureDoseSorted, self.structureDistancesSortedIndices)
        if savePredDose:
            self.predDoseVector = np.zeros(self.data.nVox)
            for s in range(self.data.nStructures):
                # with prioritization
                sReal = self.data.structureNamesInverse[self.data.weightPriorityDict[s]]
                self.predDoseVector[self.data.structureVoxels[sReal]] = self.structureDoseSorted[sReal]
                # no prioritization
                # self.predDoseVector[self.data.structureVoxels[s]] += self.structureDoseSorted[s]

    def buildMask(self):
        self.distanceMask = np.ones(self.data.nVoxFull)
        for pname in self.data.PTVNames:
            self.distanceMask[self.data.structureVoxelsFull[self.data.structureNamesInverse[pname]]] = 0
        self.distanceMask = self.distanceMask.reshape(self.data.voxDim, order='F')



    def calcAndSortDistances(self):
        self.distancesGrid = spn.distance_transform_edt(self.distanceMask, sampling=self.voxelSampling[:])
        self.distancesFull = self.distancesGrid.flatten(order='F').copy()
        self.distances = self.distancesFull[self.data.voxelIndicesInFull]
