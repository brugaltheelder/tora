from helperFunctions import *

import pickle
import numpy as np
import matplotlib.pyplot as plt
import mkl











def plotSlice(gamma,coords,doseRef, dosePred, slice):
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










f = open('predDoseOut.pkl','r')
dDict = pickle.load(f)
f.close()





assert(isinstance(dDict,dict))

print dDict['pred'].shape
print dDict.keys()
print dDict['xcoord'].shape


resolution = (0.300338,0.300338,0.25)
signed = True

from npgamma import calc_gamma

coords = (dDict['xcoord'],dDict['ycoord'], dDict['zcoord'])

print dDict['xcoord'],dDict['ycoord'],dDict['zcoord']

print np.sum(np.absolute(dDict['ref']-dDict['pred']))
print dDict['ref'].mean()





distance_threshold = 3.
distance_step_size = distance_threshold / 10.
dose_threshold = 0.03 * np.max(dDict['ref'])
lower_dose_cutoff = np.max(dDict['pred']) * 0.2
maximum_test_distance = distance_threshold * 2.


gamma = calc_gamma(coords, dDict['ref'],coords,dDict['pred'],distance_threshold,dose_threshold,lower_dose_cutoff=lower_dose_cutoff,distance_step_size=distance_step_size,maximum_test_distance=maximum_test_distance )

print gamma.shape


valid_gamma = gamma[~np.isinf(gamma)]

print valid_gamma.shape

# for i in range(0,len(dDict['xcoord']),35):
#     plt.imshow(gamma[:,i,:])
#     plt.show()


plt.hist(valid_gamma,50)
plt.show()

plotSlice(gamma,coords,dDict['ref'],dDict['pred'],60)