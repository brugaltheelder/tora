from model import *
from data import *
from predict import *
from datainputs import *
import itertools

d['modelType'] = 'fmo'

# refdosestringtag = 'medBpen'
# refdosestringtag = 'highBpen'
# refdosestringtag = 'neutral'
refdosestringtag = 'medRpen'
# refdosestringtag = 'highRpen'

refdosestring = 'refOut_' + refdosestringtag + '.mat'



# load reference data
a = sio.loadmat(d['workingDir'] + refdosestring)
doseRef = a['doseRef'].flatten()
print doseRef.shape

refdosestringtag = 'medRpenRbetter100_10'

# set reference dose
doseRefBase = doseRef.copy()

# build loop list
# dosescalar, penaltyBool, penaltyOAR, penaltyPTV
doseScalars = [1, 0.8, 1.2, 0.5, 1.5]
penalty = [(True, 1., 10.), (True, 1., 1.), (True, 0.1, 0.1), (True, 10., 10.), (False, 0., 0.)]
# penalty = [(True, 10., 1.), (True, 10., 0.1)]
iterations = [5, 10, 15]


# doseScalars = [0.8]
# penalty=[(True,1.,1.)]
# iterations = [10]

combinations = list(itertools.product(doseScalars, penalty, iterations))

for runElement in combinations:
    dScalar = runElement[0]
    penaltyBool = runElement[1][0]
    penaltyOAR = runElement[1][1]
    penaltyPTV = runElement[1][2]
    numIter = runElement[2]
    outtag = refdosestringtag + '_' + str(int(10 * dScalar)).zfill(2) + '_' + str(penaltyBool) + '_' + str(int(10 * penaltyOAR)).zfill(3) + '_' + str(int(10 * penaltyPTV)).zfill(3) + '_' + str(numIter).zfill(2)
    print outtag
    print dScalar, penaltyBool, penaltyOAR, penaltyPTV, refdosestringtag, numIter

    dat = data_fmo(d)
    doseRef = doseRefBase.copy()
    for s in range(dat.nStructures):
        # if dat.structureNames[s] not in dat.PTVNames:
        if dat.structureNames[s] in ['Rectum']:
            doseRef[dat.structureVoxels[s]] *= dScalar

    pred = predictDoseDistance(dat, doseRef)
    pred.genStructureWeightingFromArea()
    pred.predictDose()

    mod = model_fmo(dat)

    mod.data.basePenalty = penaltyBool
    mod.data.basePenaltyWeightOAR = penaltyOAR
    mod.data.basePenaltyWeightPTV = penaltyPTV
    mod.data.buildBasePenaltyVectors()

    mod.doseindex = outtag + '_00'

    plotScaled = False
    if dScalar != 1.:
        plotScaled = True

    mod.finaldoseDict['baseCase'] = doseRefBase
    mod.finaldoseDict['scaled'] = doseRef

    mod.solve()
    plotSpecific = [mod.doseindex] + ['baseCase']
    if plotScaled:
        plotSpecific += ['scaled']
    mod.pltAllDVH(saveName=mod.doseindex, plotSpecific=plotSpecific)

    n = float(numIter)
    exponentList = [0] + [1. + 1. * (i) / n for i in range(1, int(n)+1)]
    for i in range(1, int(n) + 1):
        dose = mod.finaldose.copy()
        mod.doseindex = outtag + '_' + str(i).zfill(2)
        exponent = exponentList[i]
        pred.genStructureWeightingFromArea(sourceDose=dose, ptvOverdoseExponent=exponent, oarExponent=exponent, ptvUnderdoseExponent=exponent,scaleDict={'Rectum':100.,'PTV_68':10.,'PTV_56':10.})

        newsortingdose = []
        newsortingindices = []

        pred.genDosePerStruct(dose, newsortingdose)
        pred.genSortedIndicesPerStruct(newsortingdose, newsortingindices, sortDirection=-1)
        pred.updateThreshAndWeights(pred.structureDoseSorted, newsortingindices, updateTargetDose=True, targetScalingFactor=1. - (1. * (i) / n), dualThresh=True)

        mod.solve()
        plotSpecific = [mod.doseindex] + ['baseCase']
        if plotScaled:
            plotSpecific += ['scaled']
        mod.pltAllDVH(saveName=mod.doseindex, plotSpecific=plotSpecific)

    mod.outputDose()



        # set doseindex, predictor,  then initial run then loop


        # outputs
