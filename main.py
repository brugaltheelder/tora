from model import *
from data import *
from predict import *
from datainputs import *

d['modelType'] = 'fmo'
datfmo = data_fmo(d)
mod_fmo = model_fmo(datfmo)
mod_fmo.doseindex = 'testFMO'
mod_fmo.solve()
mod_fmo.outputDose()
mod_fmo.plotDVH()
mod_fmo.pltAllDVH()

doseRef = mod_fmo.finaldose.copy()
dose = None

for o in range(datfmo.nStructures):
    if datfmo.structureNames[o] not in datfmo.PTVNames:
        # doseRef[datfmo.structureVoxels[o]] *= 0.9
        pass



d['modelType'] = 'fmo'
d['comparisonDose'] = doseRef
d['vmatFlag'] = False
d['aperLimitFlag'] = 60
dat = data_fmo(d)
predictor2 = predictDoseDistance(dat, doseRef)
predictor2.genStructureWeightingFromArea()
predictor2.predictDose()

mod = model_fmo(dat)
mod.doseindex = 'fmo'
mod.solve()
mod.outputDose()
mod.plotDVH()
mod.pltAllDVH()

n = 5.0
exponentList = [1. + 1. * (i + 1) / n for i in range(int(n))]
print exponentList
for i in range(int(n)):
    dose = mod.finaldose.copy()
    mod.doseindex = str(i)
    plotSpecific = ['original'] + [str(i)]
    exponent = exponentList[i]
    predictor2.genStructureWeightingFromArea(ptvOverdoseExponent=exponent, ptvUnderdoseExponent=exponent,
                                             oarExponent=exponent)
    newsortingdose = []
    newsortingindices = []
    predictor2.genDosePerStruct(dose, newsortingdose)
    predictor2.genSortedIndicesPerStruct(newsortingdose, newsortingindices, sortDirection=-1)
    predictor2.updateThreshAndWeights(predictor2.structureDoseSorted, newsortingindices, updateTargetDose=True,
                                      targetScalingFactor=1. - (1. * (i + 1) / n))
    mod.solve()
    mod.plotDVH(saveName=str(i))
    mod.pltAllDVH(saveName=str(i), plotSpecific=plotSpecific)
    #mod.pltAllDVH(saveName=str(i) + 'full', plotSpecific=plotSpecific, plotFull=True)


dose = mod.finaldose.copy()

dVMAT = d.copy()

dVMAT['nbeams']=180
dVMAT['beams'] = [int(360./dVMAT['nbeams']*i) for i in range(dVMAT['nbeams'])]
dVMAT['modelType'] = 'vmat'
dVMAT['comparisonDose'] = doseRef
dVMAT['vmatFlag'] = True
dVMAT['aperLimitFlag'] = None

dat = data_dao(dVMAT)
predictor = predictDoseDistance(dat, doseRef)
predictor.genStructureWeightingFromArea()
predictor.predictDose()

mod = model_dao(dat)
mod.doseindex = 'vmat'
mod.solve()
mod.outputDose()
mod.plotDVH()
mod.pltAllDVH()

dose = mod.finaldose.copy()
#dose = doseRef.copy()
n = 5.0
exponentList = [1. + 1. * (i + 1) / n for i in range(int(n))]
print exponentList
for i in range(int(n)):
    if i > 0:
        dose = mod.finaldose.copy()

    mod.doseindex = 'vmat' + str(i)
    plotSpecific = ['original'] + [mod.doseindex]
    exponent = exponentList[i]
    predictor.genStructureWeightingFromArea(ptvOverdoseExponent=exponent, ptvUnderdoseExponent=exponent,
                                            oarExponent=exponent)
    newsortingdose = []
    newsortingindices2 = []

    predictor.genDosePerStruct(dose, newsortingdose)
    predictor.genSortedIndicesPerStruct(newsortingdose, newsortingindices2, sortDirection=-1)
    predictor.updateThreshAndWeights(predictor.structureDoseSorted, newsortingindices2, updateTargetDose=True,
                                     targetScalingFactor=1. - (1. * (i + 1) / n))
    mod.solve()
    mod.plotDVH(saveName=mod.doseindex)
    mod.pltAllDVH(saveName=mod.doseindex, plotSpecific=plotSpecific)
