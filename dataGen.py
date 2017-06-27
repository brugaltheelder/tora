import sys
import os


# This is a section for the prostate case
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
voxSampling =[0.300338,0.300338,0.25]
zScale = 10.
vmatColStep = 1.0
vmatRowStep = 1.0
leafDistPerDeg = 2.0
ptvnames = ['PTV_56','PTV_68']

nbeams = 10
beams = [int(360./nbeams*i) for i in range(nbeams)]
beamsC = [0]*len(beams)
dataDir = 'Prostate/'
workingDir = dataDir + 'CASE/beams_'+str(nbeams)+'/'
vmatFlag = False
aperLimitFlag = None
useGPUFlag = False
modelType = 'fmo'
basePenaltyWeightOAR = 0.000


#
# # This is a section for the liver case
# weightUnderOverDict ={
# 'Celiac': (1.,10.),
# 'CTV': (1.,1.),
# 'DoseFalloff': (1.,10.),
# 'duodenum': (1.,10.),
# 'entrance': (1.,10.),
# 'GTV': (1.,1.),
# 'Heart': (1.,50.),
# 'KidneyL': (1.,10.),
# 'KidneyR': (1.,10.),
# 'LargeBowel': (1.,10.),
# 'Liver': (1.,20.),
# 'PTV': (500.,200.),
# 'Skin': (1.,10.),
# 'SmallBowel': (1.,10.),
# 'SMASMV': (1.,10.),
# 'SpinalCord': (1.,20.),
# 'Stomach': (1.,10.),
# }


# # weightPriorityDict = {
# # 2: 'Celiac',
# # 14: 'CTV',
# # 5: 'DoseFalloff',
# # 7: 'duodenum',
# # 0: 'entrance',
# # 15: 'GTV',
# # 12: 'Heart',
# # 3: 'KidneyL',
# # 4: 'KidneyR',
# # 9: 'LargeBowel',
# # 13: 'Liver',
# # 16: 'PTV',
# # 1: 'Skin',
# # 8: 'SmallBowel',
# # 6: 'SMASMV',
# # 10: 'SpinalCord',
# # 11: 'Stomach',
# # }

# weightPriorityDict = {
# 3: 'Celiac',
# 14: 'CTV',
# 0: 'DoseFalloff',
# 7: 'duodenum',
# 1: 'entrance',
# 15: 'GTV',
# 12: 'Heart',
# 4: 'KidneyL',
# 5: 'KidneyR',
# 9: 'LargeBowel',
# 13: 'Liver',
# 16: 'PTV',
# 2: 'Skin',
# 8: 'SmallBowel',
# 6: 'SMASMV',
# 10: 'SpinalCord',
# 11: 'Stomach',
# }

# threshDict = {
# 'Celiac': 0.,
# 'CTV': 0.,
# 'DoseFalloff': 0.,
# 'duodenum': 0.,
# 'entrance': 0.,
# 'GTV': 0.,
# 'Heart': 0.,
# 'KidneyL': 0.,
# 'KidneyR': 0.,
# 'LargeBowel': 0.,
# 'Liver': 0.,
# 'PTV': 54.,
# 'Skin': 0.,
# 'SmallBowel': 0.,
# 'SMASMV': 0.,
# 'SpinalCord': 0.,
# 'Stomach': 0.,
# }

# voxDimTup = (217, 217, 168)
# voxSampling =[0.300338,0.300338,0.25]
# zScale = 5.
# vmatColStep = 1.0
# vmatRowStep = 1.0
# leafDistPerDeg = 2.0
# ptvnames = ['PTV']
# beams = [32,74,106,148,212,254,286,328]
# beamsC = [0]*len(beams)
# # beams = [18,46,99,144,210,244,286,330]
# # beamsC = [32,13,15,32,-59,17,0,59]


# nbeams = len(beams)
# dataDir = 'Liver/'
# workingDir = dataDir + 'CASE/beams_'+str(nbeams)+'/'
# vmatFlag = False
# aperLimitFlag = None
# useGPUFlag = False
# modelType = 'fmo'
# basePenaltyWeightOAR = 0.000








d = {'modelType': modelType, 'weightUnderOverDict': weightUnderOverDict, 'weightPriorityDict': weightPriorityDict,
     'threshDict': threshDict, 'voxDim': voxDimTup, 'leafDistPerDeg': leafDistPerDeg, 'vmatColStep': vmatColStep,'voxSampling':voxSampling,'zSamplingScale':zScale,
     'vmatRowStep': vmatRowStep, 'nbeams':nbeams, 'dataDir':dataDir, 'workingDir':workingDir, 'beamsC':beamsC, 'beams':beams, 'ptvnames':ptvnames,
     'vmatFlag':vmatFlag, 'aperLimitFlag':aperLimitFlag, 'threshOverride':None, 'overWOverride':None, 'underWOverride':None,'comparisonDose':None,'useGPU':useGPUFlag, 'dualThresh':True,'basePenaltyWeightOAR':basePenaltyWeightOAR}



