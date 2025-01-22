import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib

matplotlib.use("Agg")
c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)


# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
#print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with
matplotlib.rcParams['figure.figsize'] = (15,10)
matplotlib.rcParams['lines.markersize'] = 1
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['legend.loc'] = "upper right"

pathToData = '/data/a.saricaoglu/Runs/Kroupa'
runs= [x[0] for x in os.walk(pathToData) if "COMPAS" in x[0]]
data_outputs = []
for run in runs[:2]:
    out = [f for f in os.listdir(run) if ".h5" in f]
    try:
        data = h5.File(run + "/" + out[0])
        print("Reading file: ", run + "/" + out[0])
        data_outputs.append(data)
    except:
        continue

Data = data_outputs[0]
print(list(Data.keys()))

SPs = Data ['BSE_System_Parameters']
CEs =  Data['BSE_Common_Envelopes']
MTs = Data['BSE_RLOF']
DCs = Data['BSE_Double_Compact_Objects']
RDs = Data['Run_Details']
STELLARTYPEPSTCE2 = CEs["Stellar_Type(1)"][()]
# print(RDs.keys())
# do = RDs['detailed-output']
# print(do)
mtevs =  MTs['MT_Event_Counter'][()]
firstmts = [mtevs == 1]
maskchangeMT = np.zeros(len(MTs['Stellar_Type(1)<MT'][()]), dtype='bool')
maskchangeMT2 = np.zeros(len(MTs['Stellar_Type(1)<MT'][()]), dtype='bool')

# Chooses first MT events.
maskfirstMT = (mtevs == 1)
# Chooses last MT events
masklastMT = np.zeros(len(MTs['SEED'][()]), dtype='bool')
last = 0
for i in range(0, len(CEs['SEED'][()][:20])):
    if (MTs['SEED'][()][i] == last):
        print('prev. masklastMT before switch: ', masklastMT[i-1])
        masklastMT[i-1] = False
        masklastMT[i] = True
        print('same with previous, switch')
        print('prev. seed: ', MTs['SEED'][()][i-1])
        print('current seed:', MTs['SEED'][()][i])
        print('prev. masklastMT after switch: ', masklastMT[i-1])

    else:
        masklastMT[i] = True
        last = MTs['SEED'][()][i]
        print('different from previous, no switch')
        print('prev. seed: ', MTs['SEED'][()][i-1])
        print('current seed:', MTs['SEED'][()][i])

print(MTs['SEED'][()][:20])
print(masklastMT[:20])

for i in range(0, len(MTs['Stellar_Type(1)<MT'][()][:20])):
    if (MTs['Stellar_Type(1)<MT'][()][i] == last):
        print('prev. maskchangeMT before switch: ', maskchangeMT[i-1])
        maskchangeMT[i] = False
        print('same with previous, assign false')
        print('prev. type: ', MTs['Stellar_Type(1)<MT'][()][i-1])
        print('current type:', MTs['Stellar_Type(1)<MT'][()][i])
        print('prev. maskchangeMT after switch: ', maskchangeMT[i-1])

    else:
        maskchangeMT[i] = True
        last = MTs['Stellar_Type(1)<MT'][()][i]
        print('different from previous, assign true')
        print('prev. seed: ', MTs['Stellar_Type(1)<MT'][()][i-1])
        print('current seed:', MTs['Stellar_Type(1)<MT'][()][i])

for i in range(0, len(MTs['Stellar_Type(1)<MT'][()][:20])):
    if (MTs['Stellar_Type(1)<MT'][()][i] == MTs['Stellar_Type(1)>MT'][()][i]):
        print('prev. maskchangeMT before switch: ', maskchangeMT2[i])
        maskchangeMT2[i] = False
        print('same with previous, assign false')
        print('prev. type: ', MTs['Stellar_Type(1)<MT'][()][i])
        print('current type:', MTs['Stellar_Type(1)>MT'][()][i])


    else:
        maskchangeMT2[i] = True
        print('different from previous, assign true')
        print('prev. type: ', MTs['Stellar_Type(1)<MT'][()][i])
        print('current type:', MTs['Stellar_Type(1)>MT'][()][i])

print(MTs['Stellar_Type(1)<MT'][()][:20])
print(MTs['Stellar_Type(1)>MT'][()][:20])
print(MTs['SEED'][()][:20])
print(maskchangeMT[:20])
print(maskchangeMT2[:20])
print(mtevs[:20])
# print(len(MTs['SEED']))
# print(len(np.unique(MTs['SEED'])))

# print(np.sum(firstmts))
print(SPs.keys())
print(CEs.keys())
# print(MTs.keys())

# print(len( np.unique(SPs['SEED'][()])))
# print(len(SPs['SEED'][()]))
# print(STELLARTYPEPSTCE2)
# print(str(np.sum([STELLARTYPEPSTCE2 == 13])))
m = []
mt = MTs['SEED'][()]
m.append(mt)
print(len(m))

a = np.asarray([1,2,3,1,2,3,4])
b = np.asarray([2,3,4,0,1,2,3])
c = [1,2,3,1,2,3,4]
print(np.maximum(a,b))
print(a <= 3)