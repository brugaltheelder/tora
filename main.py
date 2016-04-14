
from model import *
from data import *
from predict import *
from datainputs import *







d['modelType'] = 'fmo'
datfmo = data_fmo(d)
mod_fmo = model_fmo(datfmo)
mod_fmo.doseindex='testFMO'
mod_fmo.solve()
mod_fmo.outputDose()
mod_fmo.plotDVH()
mod_fmo.pltAllDVH()

dose = mod_fmo.finaldose

#
# d['modelType'] = 'dao'
# d['comparisonDose'] = dose
# d['vmatFlag'] = False
# d['aperLimitFlag']=60
# dat = data_dao(d)
# predictor = predictDoseDistance(dat, dose)
# predictor.genStructureWeightingFromArea()
# predictor.predictDose()
#
#
# mod = model_dao(dat)
# mod.doseindex='testVMAT'
# mod.solve()
# mod.outputDose()
# mod.plotDVH()
# mod.pltAllDVH()
#
#
# exponent = 1.0
# for i in range(3):
#     dose = mod.finaldose
#     mod.doseindex = str(i)
#     plotSpecific = ['original','testVMAT']+[str(i)]
#     exponent = 1.4*exponent
#     predictor.genStructureWeightingFromArea(ptvOverdoseExponent=exponent,ptvUnderdoseExponent=exponent,oarExponent=exponent)
#     predictor.genDosePerStruct(dose, predictor.structureDistances)
#     predictor.genSortedIndicesPerStruct(predictor.structureDistances,predictor.structureDistancesSortedIndices)
#     predictor.updateThreshAndWeights(predictor.structureDoseSorted,predictor.structureDistancesSortedIndices,updateTargetDose=True,targetScalingFactor=0.5)
#     mod.solve()
#     mod.plotDVH(saveName=str(i))
#     mod.pltAllDVH(saveName=str(i),plotSpecific=plotSpecific)


d['modelType'] = 'fmo'
d['comparisonDose'] = dose
d['vmatFlag'] = False
d['aperLimitFlag']=60
dat = data_fmo(d)
predictor2 = predictDoseDistance(dat, dose)
predictor2.genStructureWeightingFromArea()
predictor2.predictDose()


mod = model_fmo(dat)
mod.doseindex='fmo'
mod.solve()
mod.outputDose()
mod.plotDVH()
mod.pltAllDVH()



n = 5.0
exponentList = [1. + 1.*(i+1)/n for i in range(int(n))]
print exponentList
for i in range(int(n)):
    dose = mod.finaldose
    mod.doseindex = str(i)
    plotSpecific = ['original']+[str(i)]
    exponent = exponentList[i]
    predictor2.genStructureWeightingFromArea(ptvOverdoseExponent=exponent,ptvUnderdoseExponent=exponent,oarExponent=exponent)
    newsortingdose = []
    newsortingindices = []
    predictor2.genDosePerStruct(dose, newsortingdose)
    predictor2.genSortedIndicesPerStruct(newsortingdose,newsortingindices,sortDirection=-1)
    predictor2.updateThreshAndWeights(predictor2.structureDoseSorted,newsortingindices,updateTargetDose=True,targetScalingFactor=1.-(1.*(i+1)/n))
    mod.solve()
    mod.plotDVH(saveName=str(i))
    mod.pltAllDVH(saveName=str(i),plotSpecific=plotSpecific)