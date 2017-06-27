__author__ = 'troy'


import matplotlib.pyplot as plt
import numpy as np


def quadpen(dose, alpha, thresh, over=True):
    diff = dose-thresh
    if over and diff>0:
        return alpha * (diff**2)
    if not over and diff < 0:
        return alpha * (diff**2)
    return 0

def getThreshold(thresh, target, gamma):
    return gamma * thresh + (1-gamma) * target


thickness = 2


minDose = 60
maxDose = 80
doseRes = 0.01
primaryAlpha = 50.
secondaryAlpha = 1
thresh = 70.
target = 67.
gamma = 1.0

x = np.arange(minDose,maxDose,doseRes)
primaryPenalty = [0]*len(x)
secondaryPenalty = [0]*len(x)

for d in range(len(x)):
    # primary
    primthresh = getThreshold(thresh,target,gamma)
    if x[d]>max(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - max(primthresh,thresh))**2
    elif x[d]<min(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - min(primthresh,thresh))**2

    secondaryPenalty[d] = secondaryAlpha * (x[d] - thresh)**2


plt.plot(np.array(x),np.array(primaryPenalty),label = "Primary Penalty" ,linewidth=thickness)
plt.plot(np.array(x),np.array(secondaryPenalty),label = "Secondary Penalty" ,linewidth=thickness)
plt.plot((target,target),(0,5000),linewidth=thickness)
plt.plot((thresh,thresh),(0,5000),linewidth=thickness)


plt.text(70.1,4500,r'Prescription',fontsize=20)
plt.text(64.5,4500,r'Target',fontsize=20)
plt.text(74.5,3000,r'Primary',fontsize=20)
plt.text(75,200,r'Secondary',fontsize=20)

plt.xlabel("Dose")
plt.ylabel("Objective function penalty")
# plt.title("PTV Objective")
# plt.title("Initial TORA PTV Objective")
plt.savefig('objPen1.eps',dpi = 400, bbox_inches='tight')
plt.clf()



gamma = 0.0
x = np.arange(minDose,maxDose,doseRes)
primaryPenalty = [0]*len(x)
secondaryPenalty = [0]*len(x)

for d in range(len(x)):
    # primary
    primthresh = getThreshold(thresh,target,gamma)
    if x[d]>max(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - max(primthresh,thresh))**2
    elif x[d]<min(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - min(primthresh,thresh))**2

    secondaryPenalty[d] = secondaryAlpha * (x[d] - thresh)**2

plt.plot(np.array(x),np.array(primaryPenalty),label = "Primary Penalty" ,linewidth=thickness)
plt.plot(np.array(x),np.array(secondaryPenalty),label = "Secondary Penalty" ,linewidth=thickness)
plt.plot((target,target),(0,5000),linewidth=thickness)
plt.plot((thresh,thresh),(0,5000),linewidth=thickness)
plt.text(70.1,4500,r'Prescription',fontsize=20)
plt.text(64.5,4500,r'Target',fontsize=20)
plt.text(74.5,3000,r'Primary',fontsize=20)
plt.text(75,200,r'Secondary',fontsize=20)

plt.xlabel("Dose")
plt.ylabel("Objective function penalty")
# plt.title("Final TORA PTV Objective")
plt.savefig('objPen2.eps',dpi = 400, bbox_inches='tight')
plt.clf()

minDose = 0
maxDose = 40
doseRes = 0.01
primaryAlpha = 50.
secondaryAlpha = 0.5
thresh = 0.
target = 20.
gamma = 0.0

x = np.arange(minDose,maxDose,doseRes)
primaryPenalty = [0]*len(x)
secondaryPenalty = [0]*len(x)

for d in range(len(x)):
    # primary
    primthresh = getThreshold(thresh,target,gamma)
    if x[d]>max(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - max(primthresh,thresh))**2
    elif x[d]<min(primthresh,thresh):
        primaryPenalty[d] = primaryAlpha * (x[d] - min(primthresh,thresh))**2

    secondaryPenalty[d] = secondaryAlpha * (x[d] - thresh)**2

plt.plot(np.array(x),np.array(primaryPenalty),label = "Primary Penalty" ,linewidth=thickness)
plt.plot(np.array(x),np.array(secondaryPenalty),label = "Secondary Penalty" ,linewidth=thickness)
plt.plot((target,target),(0,20000),linewidth=thickness)
plt.plot((thresh+0.1,thresh+0.1),(0,20000),linewidth=thickness)
plt.text(0.1,10000,r'Prescription',fontsize=20)
plt.text(20.1,10000,r'Target',fontsize=20)
plt.text(29.5,12000,r'Primary',fontsize=20)
plt.text(30,1000,r'Secondary',fontsize=20)


plt.xlabel("Dose")
plt.ylabel("Objective function penalty")
# plt.title("OAR Objective")
# plt.title("TORA OAR Objective")

# plt.show()
plt.savefig('oarPen1.eps',dpi = 400, bbox_inches='tight')
plt.clf()