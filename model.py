try:
    import mkl

    nthreads = 4
    mkl.set_num_threads(4)
    have_mkl = True
    print("Running with MKL Acceleration with", nthreads, "cores")
except ImportError:
    have_mkl = False
    print("Running with normal backends")

import scipy.optimize as spo
import scipy.sparse as sps
from time import time
import matplotlib.pyplot as plt
import numpy as np
import os
from data import *

import numexpr as ne


class model(object):
    def __init__(self, dat, showPlots=False):
        self.data = None
        if not isinstance(dat, data):
            print 'data type wrong'
            exit()
        self.finaldoseDict = {}
        self.finalobjDict = {}
        self.doseindex = '0'
        self.showPlots = showPlots
        self.currentDose = np.zeros(dat.nVox)
        self.finaldose = None
        self.obj = 0

    def solve(self):
        pass

    def quadSolver(self, calcObjGrad, numY, x0, UB=None, display=False, ftol=1e-5, gtol=1e-5):
        res = spo.minimize(calcObjGrad, x0=x0, method='L-BFGS-B', jac=True, bounds=[(0, UB) for i in
                                                                                    xrange(
                                                                                        numY)],
                           options={'ftol': ftol, 'gtol':gtol, 'disp': display})
        return res['x'], res['fun']

    def calcObjGrad(self, yVec):
        print 'write the modality specific one for this'

    def outputDose(self):
        fullDoseOut = np.zeros(self.data.nVoxFull)
        fullDoseOut[self.data.voxelIndicesInFull] = self.finaldose.copy()
        fullDoseOut3D = fullDoseOut.reshape(self.data.voxDim, order='F')
        outputDict = {'obj': self.obj, 'dose': self.finaldose, 'fullResDose': fullDoseOut,
                      'fullResDose3d': fullDoseOut3D}

        outputDict.update(self.outputDoseModality())

        sio.savemat(self.data.workingDir + 'doseout_' + self.data.modelType + '_' + self.doseindex + '.mat',
                    outputDict)

    def outputDoseModality(self):
        pass

    def plotDVH(self, plotFull=False, saveName='', showPlots=False, saveDVH=True):
        dose = self.finaldose.copy()
        plt.clf()
        for s in range(self.data.nStructures):
            if plotFull:
                sVoxels = self.data.structureVoxelsOverlap[s]
            else:
                if not self.data.structureVoxels[s].size > 0:
                    continue
                sVoxels = self.data.structureVoxels[s]
            hist, bins = np.histogram(dose[sVoxels], bins=100)
            dvh = 1. - np.cumsum(hist) / float(sVoxels.shape[0])
            plt.plot(bins[:-1], dvh, label=self.data.structureNames[s], linewidth=2)
        lgd = plt.legend(fancybox=True, framealpha=0.5, bbox_to_anchor=(1.05, 1), loc=2)
        plt.title('DVH for doseindex:' + self.doseindex)

        plt.xlabel('Dose')
        plt.ylabel('Fractional Volume')
        if len(saveName) > 1 or saveDVH:
            saveTag = self.getSaveTag()
            plt.savefig(self.data.workingDir + 'dvh_' + saveName + '_' + saveTag + '.png',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
        if showPlots:
            plt.show()

    def getSaveTag(self):
        return 'WRITESAVETAGFUNCTION'


    def getObjUpdateZhat(self,zhat,uT,oT,uW,oW):
        oDose, uDose = np.array(self.currentDose - oT).clip(0), np.array(self.currentDose - uT).clip(-1e10,0)
        obj = float(oW.dot(oDose**2)+uW.dot(uDose**2))
        zhat += np.multiply(oW,oDose) + np.multiply(uW,uDose)
        return obj


    def pltAllDVH(self, saveName='', plotSpecific=[], saveDVH=True, plotFull=False):
        rainbow = ['r', 'c', 'darkblue', 'maroon', 'black', 'gray', 'g', 'peru', 'yellow', 'salmon', 'cadetblue']
        styles = ['solid', 'dashed', 'dotted', 'dashdot']
        plt.clf()
        title = ''
        count, lgd = 0, None

        for doseIndex, DoseDist in sorted(self.finaldoseDict.iteritems()):
            if len(plotSpecific) > 0 and doseIndex not in set(plotSpecific):
                continue

            dose = DoseDist[:]
            title = title + 'DI: ' + doseIndex + ' - ' + styles[count % len(styles)] + ' | '

            for s in range(self.data.nStructures):

                if self.data.structureVoxels[s].shape[0] > 0:
                    doseindices = None
                    if plotFull:
                        doseindices = self.data.structureVoxelsOverlap[s]
                    else:
                        doseindices = self.data.structureVoxels[s]

                    hist, bins = np.histogram(dose[doseindices], bins=100)
                    dvh = 1. - np.cumsum(hist) / float(self.data.structureVoxels[s].shape[0])
                    plt.plot(bins[:-1], dvh, label=self.data.structureNames[s], color=rainbow[s % len(rainbow)],
                             linestyle=styles[count % len(styles)], linewidth=2)
            if count == 0:
                lgd = plt.legend(fancybox=True, framealpha=0.5, bbox_to_anchor=(1.05, 1), loc=2)
            count += 1
        plt.title(title)
        maxDose = 10
        for s in range(self.data.nStructures):
            if self.data.structureNames[s] in self.data.PTVNames:
                if self.data.threshDict[self.data.structureNames[s]]*1.2>maxDose:
                    maxDose = self.data.threshDict[self.data.structureNames[s]]*1.2
        plt.axis([0,maxDose,0,1])
        plt.xlabel('Dose')
        plt.ylabel('Fractional Volume')

        if len(saveName) > 1 or saveDVH:
            saveTag = self.getSaveTag()
            plt.savefig(self.data.workingDir + 'dvh_all_' + saveName + '_' + saveTag + '.png',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
        if self.showPlots:
            plt.show()


class model_fmo(model):
    def __init__(self, dat):
        super(model_fmo, self).__init__(dat)
        if isinstance(dat, data_fmo):
            self.data = dat
            print 'data type correct for specific type'
        else:
            print 'data type wrong'
            exit()
        self.x0pB = [np.zeros(self.data.nBPB[b]) for b in range(self.data.nBeams)]
        self.x0 = np.zeros(self.data.nBeamlets)

        if self.data.comparisonDose is not None:
            self.finaldoseDict['original'] = self.data.comparisonDose

    def solve(self,ftol=1e-5, gtol=1e-5):

        start = time()

        self.fluence, obj = self.quadSolver(self.calcObjGrad, self.data.nBeamlets, self.x0.copy(), display=5, ftol=ftol, gtol=gtol)

        # self.res = spo.minimize(self.calcObjGrad, x0=self.x0.copy(), method='L-BFGS-B', jac=True, bounds=[(0, None) for i in
        #                                                                                                   xrange(
        #                                                                                                       self.data.nBeamlets)],
        #                         options={'ftol': 1e-4, 'disp': 5})
        print 'solved in ', time() - start, ' seconds'
        # self.fluence = self.res['x']
        self.fluencepB = [self.fluence[self.data.nBPBcum[b]:self.data.nBPBcum[b + 1]] for b in range(self.data.nBeams)]
        self.finaldose = self.currentDose.copy()
        self.finaldoseDict[self.doseindex] = self.finaldose.copy()
        self.finalobjDict[self.doseindex] = obj
        self.obj = obj

    def calcDose(self, fluence):
        self.currentDose = np.zeros(self.data.nVox)
        for b in range(self.data.nBeams):
            self.currentDose += self.data.dijList[b].dot(fluence[self.data.nBPBcum[b]:self.data.nBPBcum[b + 1]])

    def calcObjGrad(self, fluence):
        self.calcDose(fluence)

        # oDose, uDose = np.array(self.currentDose - self.data.thresh), np.array(self.currentDose - self.data.thresh)
        # oDose, uDose = np.array(self.currentDose - self.data.overThresh), np.array(self.currentDose - self.data.underThresh)
        # obj = float(
        #     self.data.overPenalty.dot(oDose.clip(0) ** 2) + self.data.underPenalty.dot(uDose.clip(-1e10, 0) ** 2))
        #
        # zhat = np.multiply(self.data.overPenalty, oDose.clip(0)) + np.multiply(self.data.underPenalty,
        #                                                                        uDose.clip(-1e10, 0))
        grad = np.zeros(self.data.nBeamlets)
        zhat = np.zeros(self.data.nVox)
        obj = 0.0
        obj += self.getObjUpdateZhat(zhat, self.data.underThresh,self.data.overThresh,self.data.underPenalty,self.data.overPenalty)

        if self.data.basePenalty:
            obj += self.getObjUpdateZhat(zhat, self.data.baseThresh,self.data.baseThresh,self.data.basePenaltyUnder,self.data.basePenaltyOver)


        for b in range(self.data.nBeams):
            grad[self.data.nBPBcum[b]:self.data.nBPBcum[b + 1]] = self.data.dijList[b].transpose().dot(zhat)

        return obj, grad

    def getSaveTag(self):
        savetag = 'fmo'
        return savetag

    def outputDoseModality(self):
        outputDict = {'obj': self.obj, 'fluence': self.fluence, 'fluencepB': self.fluencepB, 'thresh': self.data.thresh,
                      'overThresh':self.data.overThresh,'underThresh':self.data.underThresh,
                      'over': self.data.overPenalty, 'under': self.data.underPenalty}
        return outputDict


class model_dao(model):
    def __init__(self, dat):
        super(model_dao, self).__init__(dat)
        if isinstance(dat, data_dao):
            self.data = dat
            print 'data type correct for specific type'
        else:
            print 'data type wrong'
            exit()

        self.y = np.zeros(self.data.aperLimit)
        self.bIndexInY = np.zeros(self.data.aperLimit, dtype='int32')
        self.objPerIter = []
        self.yAperBeamlets = [None] * self.data.aperLimit
        self.aperColLeftPerRow = [[-1 for r in range(self.data.spacedRowMax)] for b in range(self.data.aperLimit)]
        self.aperColRightPerRow = [[-1 for r in range(self.data.spacedRowMax)] for b in range(self.data.aperLimit)]
        self.currentY = 0

        self.data.Dkj = np.zeros((self.data.nVox, self.data.aperLimit))
        self.beamUsed = np.zeros(self.data.nBeams, dtype='int32')
        if self.data.comparisonDose is not None:
            self.finaldoseDict['baseCase'] = self.data.comparisonDose

    def solveRMP(self):

        if self.currentY > 0:

            self.data.DkjRMP = sps.csc_matrix(self.data.Dkj[:, 0:self.currentY])

            self.y[0:self.currentY], obj = self.quadSolver(self.calcObjGrad, self.currentY,
                                                           self.y[0:self.currentY].copy(), UB=None)

            # self.yRMP = self.y[0:self.currentY]
            #
            # self.res = spo.minimize(self.calcObjGrad, x0=self.yRMP.copy(), method='L-BFGS-B', jac=True, bounds=[(0, None) for i in
            #                                                                                               xrange(
            #                                                                                                   self.currentY)],
            #                     options={'ftol': 1e-4, 'disp': False})
            # self.y[0:self.currentY] = self.res['x']
            self.objPerIter.append(obj)
            self.calcDose(self.y[0:self.currentY])
            self.finaldose = self.currentDose.copy()

        else:
            self.finaldose = np.zeros(self.data.nVox)
            self.objPerIter.append(self.calcObjGrad(self.y[0:self.currentY])[0])

    def calcDose(self, yVec):
        if len(yVec) > 0:
            self.currentDose = self.data.DkjRMP.dot(yVec)
        else:
            self.currentDose.fill(0.)

    def calcObjGrad(self, yVec):
        self.calcDose(yVec)
        # oDose = np.array(self.currentDose - self.data.thresh).clip(0)
        # uDose = oDose.clip(-1e10, 0)
        # oDose = oDose.clip(0)

        oDose = np.array(self.currentDose - self.data.overThresh).clip(0)
        uDose = np.array(self.currentDose - self.data.underThresh).clip(-1e10, 0)



        obj = float(self.data.overPenalty.dot(ne.evaluate('oDose ** 2')) + self.data.underPenalty.dot(
            ne.evaluate('uDose ** 2')))

        zhat = np.multiply(self.data.overPenalty, oDose) + np.multiply(self.data.underPenalty, uDose)
        if self.currentY > 0:
            grad = self.data.DkjRMP.transpose().dot(zhat)
        else:
            grad = []
        return obj, grad

    def getSaveTag(self):
        savetag = ''
        if self.data.vmatFlag:
            savetag = savetag + 'vmat_' + str(self.data.aperLimit)
        else:
            savetag = savetag + 'dao_' + str(self.data.aperLimit)
        return savetag

    def solve(self):

        beginning = time()
        # clear beams
        self.beamUsed = np.zeros(self.data.nBeams)
        self.currentY = 0
        self.currentDose = np.zeros(self.data.nVox)
        del self.objPerIter[:]

        while self.currentY < self.data.aperLimit:
            # solve RMP
            # construct y vector, construct Dij matrix

            start = time()

            self.solveRMP()

            # if self.currentY>0:
            #
            #     self.data.DkjRMP = sps.csc_matrix(self.data.Dkj[:,0:self.currentY])
            #
            #     self.y[0:self.currentY],obj = self.quadSolver(self.calcObjGrad,self.currentY,self.y[0:self.currentY].copy(),UB=None)
            #
            #     # self.yRMP = self.y[0:self.currentY]
            #     #
            #     # self.res = spo.minimize(self.calcObjGrad, x0=self.yRMP.copy(), method='L-BFGS-B', jac=True, bounds=[(0, None) for i in
            #     #                                                                                               xrange(
            #     #                                                                                                   self.currentY)],
            #     #                     options={'ftol': 1e-4, 'disp': False})
            #     # self.y[0:self.currentY] = self.res['x']
            #     self.objPerIter.append(obj)
            #     self.calcDose(self.y)
            #     self.finaldose = self.currentDose.copy()
            #
            # else:
            #     self.finaldose = np.zeros(self.data.nVox)
            #     self.objPerIter.append(self.calcObjGrad(self.y[0:self.currentY])[0])

            rpmTime = time()

            # todo do pricing problem case switch when it is ready here
            if self.data.vmatFlag:
                self.pricingProblemVMATBlocks(self.y[0:self.currentY])
            else:
                self.pricingProblemDAO(self.y[0:self.currentY])

            print 'iter:', self.currentY, 'RMP in', rpmTime - start, 's PP in ', time() - rpmTime, ' s with b', \
            self.bIndexInY[self.currentY], 'with', len(
                self.yAperBeamlets[self.currentY]), 'beamlets in the aperture and obj', self.objPerIter[-1]

            self.currentY += 1

        # final RMP solution
        self.yRMP = self.y[0:self.currentY]
        self.data.DkjRMP = sps.csc_matrix(self.data.Dkj[:, 0:self.currentY])
        self.y[0:self.currentY], obj = self.quadSolver(self.calcObjGrad, self.currentY, self.y[0:self.currentY].copy(),
                                                       UB=None)
        self.objPerIter.append(obj)
        #
        # self.res = spo.minimize(self.calcObjGrad, x0=self.yRMP.copy(), method='L-BFGS-B', jac=True, bounds=[(0, None) for i in
        #                                                                                               xrange(
        #                                                                                                   self.currentY)],
        #                     options={'ftol': 1e-4, 'disp': self.optDisplay})
        # self.y[0:self.currentY] = self.res['x']
        # self.objPerIter.append(self.res['fun'])
        self.calcDose(self.y)
        self.finaldose = self.currentDose.copy()
        self.finaldoseDict[self.doseindex] = self.currentDose.copy()
        self.finalobjDict[self.doseindex] = obj
        self.obj = obj

    def getBeamletGrad(self, zhat, b):
        return self.data.dijList[b].transpose().dot(zhat)

    def pricingProblemVMATBlocks(self, yVec):
        # best aper, best beam saving
        bestBeamIndex = -1
        bestAperBeamlets = []
        bestAperScore = 0
        bestEndpointsL = []
        bestEndpointsR = []
        # calc helpers
        self.calcDose(yVec)
        oDose, uDose = np.array(self.currentDose - self.data.overThresh), np.array(self.currentDose - self.data.underThresh)
        #oDose, uDose = np.array(self.currentDose - self.data.thresh), np.array(self.currentDose - self.data.thresh)

        zhat = np.multiply(self.data.overPenalty, oDose.clip(0)) + np.multiply(self.data.underPenalty,
                                                                               uDose.clip(-1e10, 0))

        for b in range(self.data.nBeams):

            if not self.beamUsed[b]:
                # calc grad
                grad = self.getBeamletGrad(zhat, b)
                worth = 0
                beamlets = []
                # run pricing problem for each beam
                worth, beamlets, endPointsL, endPointsR = self.pricingProblemVMATBlocksBeam(b, grad)

                if worth < bestAperScore:
                    bestAperScore = worth
                    bestAperBeamlets = beamlets[:]
                    bestBeamIndex = b
                    bestEndpointsL = endPointsL[:]
                    bestEndpointsR = endPointsR[:]

        bestBeamletBool = np.zeros(self.data.nBPB[bestBeamIndex])
        self.yAperBeamlets[self.currentY] = bestAperBeamlets[:]
        bestBeamletBool[np.ix_(np.array(bestAperBeamlets))] = 1

        self.aperColLeftPerRow[bestBeamIndex] = bestEndpointsL[:]
        self.aperColRightPerRow[bestBeamIndex] = bestEndpointsR[:]

        # UPDATE LEFT AND RIGHT COLUMNS OF THE NEW APERTURE

        self.data.Dkj[:, self.currentY] = self.data.dijList[bestBeamIndex].dot(bestBeamletBool)

        self.bIndexInY[self.currentY] = bestBeamIndex
        if self.data.vmatFlag:
            self.beamUsed[bestBeamIndex] = 1

    def pricingProblemVMATBlocksBeam(self, beamNum, gradient):

        endPointsL, endPointsR = [-1] * self.data.spacedRowMax, [-1] * self.data.spacedRowMax
        beamlets, worth = [], 0

        for r in range(len(self.data.rowStartPB[beamNum])):
            # calc best aper
            maxSoFar, maxEndingHere, lE, rE, = 0, 0, self.data.rowStartPB[beamNum][r], self.data.rowStartPB[beamNum][r]
            rowSize = self.data.rowEndPB[beamNum][r] - self.data.rowStartPB[beamNum][r]
            for i in range(self.data.rowStartPB[beamNum][r], self.data.rowEndPB[beamNum][r]):
                maxEndingHere += gradient[i]
                if maxEndingHere > 0:
                    maxEndingHere, lE, rE = 0, i + 1, i + 1
                if maxSoFar > maxEndingHere:
                    maxSoFar, rE = maxEndingHere, i + 1

            # calculate LL, LR, RL, RR
            LL, LR, RL, RR = 0, self.data.spacedRowMax, 0, self.data.spacedRowMax
            for j in range(beamNum - 1, -1, -1):
                if not self.beamUsed[j]:
                    continue
                elif self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] > -1:
                    LL = max(LL, self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] - int(
                        self.data.leafDistPerDeg * (beamNum - j) * self.data.degPerStep))
                    LR = min(LR, self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] + int(
                        self.data.leafDistPerDeg * (beamNum - j) * self.data.degPerStep))
                    RR = min(RR, self.aperColRightPerRow[j][self.data.rowIndToActual[beamNum][r]] + int(
                        self.data.leafDistPerDeg * (beamNum - j) * self.data.degPerStep))
                    RL = max(RL, self.aperColRightPerRow[j][self.data.rowIndToActual[beamNum][r]] - int(
                        self.data.leafDistPerDeg * (beamNum - j) * self.data.degPerStep))
                    break
            for j in range(beamNum + 1, self.data.nBeams):
                if not self.beamUsed[j]:
                    continue
                elif self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] > -1:
                    LL = max(LL, self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] - int(
                        self.data.leafDistPerDeg * (-beamNum + j) * self.data.degPerStep))
                    LR = min(LR, self.aperColLeftPerRow[j][self.data.rowIndToActual[beamNum][r]] + int(
                        self.data.leafDistPerDeg * (-beamNum + j) * self.data.degPerStep))
                    RR = min(RR, self.aperColRightPerRow[j][self.data.rowIndToActual[beamNum][r]] + int(
                        self.data.leafDistPerDeg * (-beamNum + j) * self.data.degPerStep))
                    RL = max(RL, self.aperColRightPerRow[j][self.data.rowIndToActual[beamNum][r]] - int(
                        self.data.leafDistPerDeg * (-beamNum + j) * self.data.degPerStep))
                    break

            if LL > LR or RR < RL:
                print beamNum
                print LL, LR, RL, RR, 'execption error'

            if lE - self.data.rowStartPB[beamNum][r] + self.data.colStartPerRow[beamNum][r] < LL:
                lE = self.data.rowStartPB[beamNum][r] - self.data.colStartPerRow[beamNum][r] + LL
            elif lE - self.data.rowStartPB[beamNum][r] + self.data.colStartPerRow[beamNum][r] > LR:
                lE = self.data.rowStartPB[beamNum][r] - self.data.colStartPerRow[beamNum][r] + LR
            if rE - self.data.rowStartPB[beamNum][r] + self.data.colStartPerRow[beamNum][r] < RL:
                rE = self.data.rowStartPB[beamNum][r] - self.data.colStartPerRow[beamNum][r] + RL
            elif rE - self.data.rowStartPB[beamNum][r] + self.data.colStartPerRow[beamNum][r] > RR:
                rE = self.data.rowStartPB[beamNum][r] - self.data.colStartPerRow[beamNum][r] + RR

            endPointsL[self.data.rowIndToActual[beamNum][r]] = lE - self.data.rowStartPB[beamNum][r] + \
                                                               self.data.colStartPerRow[beamNum][r]
            endPointsR[self.data.rowIndToActual[beamNum][r]] = rE - self.data.rowStartPB[beamNum][r] + \
                                                               self.data.colStartPerRow[beamNum][r]

            for i in range(lE, rE):
                beamlets.append(i)

            # run PP on row, return column endpoints
            worth += maxSoFar
        return worth, beamlets, endPointsL, endPointsR

    def pricingProblemDAO(self, yVec, useCurrent=True):
        # best aper, best beam saving
        bestBeamIndex = -1
        bestAperBeamlets = []
        bestAperScore = 0
        # calc helpers
        self.calcDose(yVec)
        # oDose, uDose = np.array(self.currentDose - self.data.thresh), np.array(self.currentDose - self.data.thresh)
        oDose, uDose = np.array(self.currentDose - self.data.overThresh), np.array(self.currentDose - self.data.underThresh)
        zhat = np.multiply(self.data.overPenalty, oDose.clip(0)) + np.multiply(self.data.underPenalty,
                                                                               uDose.clip(-1e10, 0))
        for b in range(self.data.nBeams):
            if not self.beamUsed[b]:
                grad = self.getBeamletGrad(zhat, b)
                worth, beamlets = self.pricingProblemBeamDAO(grad, b)
                if worth < bestAperScore:
                    bestAperScore = worth
                    bestAperBeamlets = beamlets[:]
                    bestBeamIndex = b
        bestBeamletBool = np.zeros(self.data.nBPB[bestBeamIndex])
        self.yAperBeamlets[self.currentY] = bestAperBeamlets[:]
        bestBeamletBool[np.ix_(np.array(bestAperBeamlets))] = 1

        # these worse than using a dense Dkj
        # self.Dkj.getcol(self.currentY)[:] = self.data.dijList[bestBeamIndex].dot(bestBeamletBool).reshape((self.data.nVox,1)).copy()
        # self.Dkj[:,self.currentY] = self.data.dijList[bestBeamIndex].dot(bestBeamletBool).reshape((self.data.nVox,1)).copy()
        # self.Dkj[:,self.currentY] = sps.coo_matrix(self.data.dijList[bestBeamIndex].dot(bestBeamletBool)).transpose().tocsc()
        self.data.Dkj[:, self.currentY] = self.data.dijList[bestBeamIndex].dot(bestBeamletBool)

        self.bIndexInY[self.currentY] = bestBeamIndex
        if self.data.vmatFlag:
            self.beamUsed[bestBeamIndex] = 1

    def pricingProblemBeamDAO(self, gradient, b):
        beamlets, worth = [], 0
        for l in range(len(self.data.rowStartPB[b])):
            maxSoFar, maxEndingHere, lE, rE, = 0, 0, self.data.rowStartPB[b][l], self.data.rowStartPB[b][l]
            for i in range(self.data.rowStartPB[b][l], self.data.rowEndPB[b][l]):
                maxEndingHere += gradient[i]
                if maxEndingHere > 0:
                    maxEndingHere, lE, rE = 0, i + 1, i + 1
                if maxSoFar > maxEndingHere:
                    maxSoFar, rE = maxEndingHere, i + 1
            for i in range(lE, rE):
                beamlets.append(i)
            worth += maxSoFar
        return worth, beamlets

    def outputDoseModality(self):
        outputDict = {'obj': self.obj, 'y': self.y, 'beamlets': self.yAperBeamlets, 'thresh': self.data.thresh,
                      'overThresh':self.data.overThresh,'underThresh':self.data.underThresh,
                      'over': self.data.overPenalty, 'under': self.data.underPenalty}
        return outputDict
