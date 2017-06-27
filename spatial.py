from model import *
from data import *
from predict import *
from datainputs import *
import itertools

import pickle



print 'Running spatial DVH file'


d['modelType'] = 'fmo'



# build list of methods

spatialSortingMethods = ['distanceToPTV']
TORArefinementFlag = True

# build name of run options
refdosestringtag = ''

aBoost = ['default','PTV','Bladder','Rectum']
OARscaling = [1., 0.8, 1.2]
BladderScaling = [1.,0.8,1.2]
RectumScaling = [1.,0.8,1.2]

combinations = list(itertools.product(aBoost, OARscaling, BladderScaling,RectumScaling))

iternumBoost = 0
iternumOAR = 0
iternumBlad = 0
iternumRect = 0
boost = aBoost[iternumBoost]
OARscale = OARscaling[iternumOAR]
Bladscale = BladderScaling[iternumBlad]
Rectscale = RectumScaling[iternumRect]

refdosestringtag = boost + '_' + str(int(10*OARscale)).zfill(2) + '_' + str(int(10*Bladscale)).zfill(2) + '_' + str(int(10*Rectscale)).zfill(2)

print "reference dose tag:", refdosestringtag

#generate data file
print "generating data"
dat = data_fmo(d)



#generate alphas for dose distribution
genAlphas(dat,PTVBase=100.,OARBase=1.,modifier=['Bladder'])


#find ground truth dose distribution (solve FMO)
mod_fmo_GT = model_fmo(dat)

mod_fmo_GT.doseindex = refdosestringtag + '_' + 'base'

mod_fmo_GT.solve()

#save ground truth dose distribution (maybe make folder structure)


mod_fmo_GT.outputDose()
# todo modify GT stuff here for the OAR/RECT/BLAD boosting
doseRefGT = mod_fmo_GT.finaldose.copy()



for m in range(len(spatialSortingMethods)):


    pred = None
    predDose = np.zeros(dat.nVox)
    # build predictor object
    # evaluate method spatialSortingMethods[m]

    if spatialSortingMethods[m]=='distanceToPTV':
        print 'using distance to PTV'
        pred=predictDoseDistance(dat,doseRefGT)
        pred.predictDose(savePredDose=True,updateWeights=False)
        predDose = pred.predDoseVector.copy()


    elif spatialSortingMethods[m]=='default':
        print 'put in a real value'
    else:
        print 'no predictor found'


    doseRefGT3D = np.zeros(dat.nVoxFull)
    doseRefGT3D[dat.voxelIndicesInFull] = doseRefGT
    doseRefGT3D = doseRefGT3D.reshape(dat.voxDim,order='F').copy()
    predRefGT3D = np.zeros(dat.nVoxFull)
    predRefGT3D[dat.voxelIndicesInFull] = predDose
    predRefGT3D = predRefGT3D.reshape(dat.voxDim,order='F').copy()

    print 'total dose ref:',doseRefGT.sum()
    print 'total dose pred:',predRefGT3D.sum()

    #add in coordinates here

    xc, yc, zc = pred.getSpatialGrid()


    if TORArefinementFlag:
        pass
        # run tora on inputs (build this function)
        #save outputs to output dict

    # dDict = {'ref':doseRefGT3D,'pred':predRefGT3D, 'xcoord':xc, 'ycoord':yc,'zcoord':zc }
    # pFile = open('predDoseOut.pkl','w')
    # pickle.dump(dDict,pFile)
    # pFile.close()


    # run dose evaluation

    # print evalulation, save outputs







