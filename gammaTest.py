from helperFunctions import *

import pickle
import numpy as np

f = open('predDoseOut.pkl','r')
dDict = pickle.load(f)
f.close()

distance = 4.
threshold = 2.
resolution = (0.300338,0.300338,0.25)
signed = True

gamma_map = gamma_evaluation(dDict['pred'],dDict['ref'],distance,threshold,resolution,signed)
print gamma_map