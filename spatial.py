from model import *
from data import *
from predict import *
from datainputs import *
import itertools



print 'Running spatial DVH file'


d['modelType'] = 'fmo'



# build list of methods

spatialSortingMethods = ['distance']


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
genAlphas(dat,PTVBase=1000000.,OARBase=1.,modifier=['Bladder'])



#find ground truth dose distribution (solve FMO)
mod_fmo_GT = model_fmo(dat)

mod_fmo_GT.doseindex = refdosestringtag + '_' + 'base'

mod_fmo_GT.solve()

#save ground truth dose distribution (maybe make folder structure)


mod_fmo_GT.outputDose()
# todo modify GT stuff here for the OAR/RECT/BLAD boosting
doseRefGT = mod_fmo_GT.finaldose.copy()



for m in range(len(spatialSortingMethods)):
    pass
    # build predictor object


    # read in ground truth dose


    # evaluate method spatialSortingMethods[m]


    # evaluate comparisons for the dose


    # save and output results


