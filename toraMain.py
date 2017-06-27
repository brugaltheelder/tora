# Imports

from model import *
from data import *
from predict import *
from datainputs import *
import itertools


# Helper run funcitons

def loadReferenceDose(d,refdosestring, refdosesourcename, refdosestringname, refdosebasestringname):
    a = sio.loadmat(d['workingDir']+refdosestring)
    d[refdosestringname] = a[refdosesourcename].flatten()
    d[refdosebasestringname] = d[refdosestringname].copy()

def scaleReferenceDose(dat,refDose,dScalar,scaledStructures):
    assert(isinstance(dat,data_fmo))
    for s in range(dat.nStructures):
        if dat.structureNames[s] in scaledStructures:
            doseRef[dat.structureVoxels[s]] *= dScalar







# case inputs - Prostate


refdosestringtag = 'medBpen'
# refdosestringtag = 'highBpen'
# refdosestringtag = 'neutral'
# refdosestringtag = 'medRpen'
# refdosestringtag = 'highRpen'

doseRefNameInSource = 'doseRef'

refdosestring = 'refOut_' + refdosestringtag + '.mat'


# refdosestringtag = 'paperFMO_liver'
# doseRefNameInSource = 'dose'
# refdosestring = 'doseout_fmo_' + refdosestringtag + '.mat'

loadReferenceDose(d, refdosestring, doseRefNameInSource, 'doseRef', 'doseRefBase')




# refdosestringtag = 'paperFMO_liver'
# refdosestringtag = 'sourceFMOtuning'



# doseRefNameInSource = 'dose'
# refdosestring = 'doseout_fmo_' + refdosestringtag + '.mat'

# loadReferenceDose(d, refdosestring, doseRefNameInSource, 'doseRef', 'doseRefBase')



doseScalars = [1.,1.15,0.85]
# doseScalars =[1.]
#penalty = [(True, 0.00001, 0.00001)]
penalty = [(True, 0.001, 0.001)]
# penalty = [(True, 1., 1.)]

#penalty = [(True, 1., 1.), (True, 10., 10.),(True, 0.001, 0.001)]
# iterations = [1,2,3,5]
iterations = [10]
scaledOARs = ['Rectum']
# scaledOARs = ['Liver']
dontPlot = ['SMASMV','entrance','duodenum','DoseFalloff','CTV','GTV','Celiac']
# otherTag = '_Bweightboost'
otherTag = 'EvenWeights'
ows=1.
owsn=''
oarLower, oarUpper, ptvLower, ptvUpper = 0., 1., 1., 1.
plotTitleBool = False

# generate combinations for iterating

combinations = list(itertools.product(doseScalars, penalty, iterations))


for runElement in combinations:
    # construct helper variables
    dScalar = runElement[0]
    penaltyBool = runElement[1][0]
    penaltyOAR = runElement[1][1]
    penaltyPTV = runElement[1][2]
    numIter = runElement[2]
    scaledOARtag = ''.join([i[0] for i in scaledOARs])
    outtag = refdosestringtag + '_' + str(int(10 * dScalar)).zfill(2) + '_' + str(penaltyBool) + '_' + str(int(10 * penaltyOAR)).zfill(3) + '_' + str(int(10 * penaltyPTV)).zfill(3) + '_' + str(numIter).zfill(2) + '_' + scaledOARtag + otherTag
    plotScaled = False
    if dScalar != 1.:
        plotScaled = True
    print outtag
    print dScalar, penaltyBool, penaltyOAR, penaltyPTV, refdosestringtag, numIter

    # generate data object from d
    dat = data_fmo(d)
    doseRef = d['doseRef'].copy()
    doseRefBase = d['doseRefBase'].copy()

    # scale reference dose
    scaleReferenceDose(dat, doseRef,dScalar,scaledOARs)

    # generate predictor object and initial predition
    pred = predictDoseDistance(dat, doseRef)
    #pred.genStructureWeightingFromArea()
    pred.setStructureWeightingFromValues(oarLower,oarUpper,ptvLower,ptvUpper, overrideWeightScalar=ows, overrideWeightStructureName=owsn)
    pred.predictDose()

    # generate model object from dat
    mod = model_fmo(dat)

    mod.data.basePenalty = penaltyBool
    mod.data.basePenaltyWeightOAR = penaltyOAR
    mod.data.basePenaltyWeightPTV = penaltyPTV
    mod.data.buildBasePenaltyVectors(doseScalar=doseRef)
    mod.doseindex = outtag + '_00'

    mod.finaldoseDict['baseCase'] = doseRefBase
    mod.finaldoseDict['scaled'] = doseRef

    if plotScaled:
        plotSpecific =  ['scaled'] + [mod.doseindex] + ['baseCase']
    else:
        plotSpecific = [mod.doseindex] + ['baseCase']

    # invoke initial solve of model
    mod.solve()

    mod.pltAllDVH(saveName=mod.doseindex, plotSpecific=plotSpecific,doNotPlotTheseStructures=dontPlot, titleOverride='Solid = Auto-planned DVH, Dashed = Reference DVH', plotTitle=plotTitleBool)



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
        pred.updateThreshAndWeights(pred.structureDoseSorted, newsortingindices, updateTargetDose=True, targetScalingFactor=1. - (1. * (i) / n), dualThresh=True, updateWeights=False)


        mod.solve()
        if plotScaled:
            plotSpecific = ['scaled'] + [mod.doseindex] + ['baseCase']
        else:
            plotSpecific = [mod.doseindex] + ['baseCase']
        mod.pltAllDVH(saveName=mod.doseindex, plotSpecific=plotSpecific,doNotPlotTheseStructures=dontPlot, titleOverride='Solid = Auto-planned DVH, Dashed = Reference DVH', plotTitle=plotTitleBool)

    mod.outputDose()








