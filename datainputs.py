import sys
import os
weightUnderOverDict = {'BODY': (1., 1.),
                       'Bladder': (1., 10.),
                       'Lt_femoral_head': (1., 10.),
                       'Lymph_Nodes': (1., 10.),
                       'PTV_56': (500., 200.),
                       'PTV_68': (2000., 800.),
                       'Penile_bulb': (1., 10.),
                       'Rectum': (1., 10.),
                       'Rt_femoral_head': (1., 10.),
                       'prostate_bed': (1., 10.)}
weightPriorityDict = {9: 'PTV_68',
                      6: 'Rectum',
                      3: 'Lt_femoral_head',
                      0: 'BODY',
                      7: 'prostate_bed',
                      2: 'Rt_femoral_head',
                      1: 'Lymph_Nodes',
                      8: 'PTV_56',
                      4: 'Penile_bulb',
                      5: 'Bladder'}
threshDict = {'BODY': 0.,
              'Bladder': 0.,
              'Lt_femoral_head': 0.,
              'Lymph_Nodes': 0.,
              'PTV_56': 56.,
              'PTV_68': 68.,
              'Penile_bulb': 0.,
              'Rectum': 0.,
              'Rt_femoral_head': 0.,
              'prostate_bed': 0.}

voxDimTup = (184, 184, 90)
# voxSampling =[0.300338,0.300338,0.25]
voxSampling = [0.300338, 0.300338, 2.5]
# voxSampling =[0.300338,0.300338,4]
vmatColStep = 1.0
vmatRowStep = 1.0
leafDistPerDeg = 2.0

ptvnames = ['PTV_56','PTV_68']

if len(sys.argv) > 1:
    nbeams = int(sys.argv[1])
else:
    nbeams = 10

beams = [int(360./nbeams*i) for i in range(nbeams)]


dataDir = 'Prostate/'
workingDir = dataDir + 'CASE/beams_'+str(nbeams)+'/'
if not os.path.isdir(workingDir):
    os.makedirs(workingDir)

vmatFlag = False
aperLimitFlag = None
useGPUFlag = False

d = {'modelType': 'fmo', 'weightUnderOverDict': weightUnderOverDict, 'weightPriorityDict': weightPriorityDict,
     'threshDict': threshDict, 'voxDim': voxDimTup, 'leafDistPerDeg': leafDistPerDeg, 'vmatColStep': vmatColStep,'voxSampling':voxSampling,
     'vmatRowStep': vmatRowStep, 'nbeams':nbeams, 'dataDir':dataDir, 'workingDir':workingDir, 'beams':beams, 'ptvnames':ptvnames,
     'vmatFlag':vmatFlag, 'aperLimitFlag':aperLimitFlag, 'threshOverride':None, 'overWOverride':None, 'underWOverride':None,'comparisonDose':None,'useGPU':useGPUFlag, 'dualThresh':True}


# import json
# with open('testDictOut.json','w') as f:
#     json.dump(d,f)
