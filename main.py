
from model import *
from data import *
from datainputs import *






dat = data_dao(d)
mod = model_dao(dat)
mod.doseindex='testVMAT'
mod.solve()
mod.outputDose()
mod.plotDVH()
mod.pltAllDVH()

d['modelType'] = 'fmo'
datfmo = data_fmo(d)
mod_fmo = model_fmo(datfmo)
mod_fmo.doseindex='testFMO'
mod_fmo.solve()
mod_fmo.outputDose()
mod_fmo.plotDVH()
mod_fmo.pltAllDVH()