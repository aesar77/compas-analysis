import os, sys
import numpy as np               # for handling arrays
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
import matplotlib.pyplot as plt  #for plotting

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with

pathToData = 'COMPAS_Output_2/COMPAS_Output.h5'

# This is known as an ipython magic command, and allows plots to be produced within the notebook
Data  = h5.File(pathToData)
print(list(Data.keys()))
SPs = Data['BSE_System_Parameters']
MTs = Data['BSE_RLOF']
CEs = Data['BSE_Common_Envelopes']
SNe = Data['BSE_Supernovae']
DCs = Data['BSE_Double_Compact_Objects']

print(CEs.keys())

CEpos = []
CEneg = []

seeds, events = getEventHistory(Data)
print(events)

for ii, seed in enumerate(seeds):
    # print(seed, events[ii])
    occurance = False
    for event in events[ii]:

        if event[0] == "MT" and event[-1] == True and occurance == False:
            #print(event[0], event[-1])
            CEpos.append(seed)
            occurance = True
        elif seed not in CEneg:
            CEneg.append(seed)

# print(len(CEpos))
# print(len(CEneg))
# print(len(seeds))

CEP_sys = []
CEN_sys = []
CE_sys = []

# for negseed in CEneg:
#     for mtSeedIndex in range(len(MTs['SEED'][()])):
#         mtSeed = MTs['SEED'][mtSeedIndex]
#         if mtSeed == negseed:
#             CEN_sys.append([negseed, MTs["Mass(1)<MT"][mtSeedIndex],MTs["Mass(1)>MT"][mtSeedIndex],MTs["Mass(2)<MT"][mtSeedIndex],MTs["Mass(2)>MT"][mtSeedIndex], MTs['Eccentricity<MT'][mtSeedIndex], MTs['Eccentricity>MT'][mtSeedIndex]])

# for posseed in CEpos:
#     for ceSeedIndex in range(len(CEs['SEED'][()])):
#         ceSeed = CEs['SEED'][ceSeedIndex]
#         if ceSeed == posseed:
#             CE_sys.append([posseed, CEs["Mass(1)<CE"][ceSeedIndex],CEs["Mass(1)>CE"][ceSeedIndex],CEs["Mass(2)<CE"][ceSeedIndex],CEs["Mass(2)>CE"][ceSeedIndex], CEs['Eccentricity<CE'][ceSeedIndex], CEs['Eccentricity>CE'][ceSeedIndex]])
#     for mtSeedIndex in range(len(MTs['SEED'][()])):
#         mtSeed = MTs['SEED'][mtSeedIndex]
#         if mtSeed == posseed:
#             CEP_sys.append([posseed, MTs["Mass(1)<MT"][mtSeedIndex],MTs["Mass(1)>MT"][mtSeedIndex],MTs["Mass(2)<MT"][mtSeedIndex],MTs["Mass(2)>MT"][mtSeedIndex], MTs['Eccentricity<MT'][mtSeedIndex], MTs['Eccentricity>MT'][mtSeedIndex]])

# plt.subplot(2,1,1)
# for s in CEN_sys:
#     plt.plot(s[1],s[3], marker = "o", color = "b")
# for s in CEP_sys:
#     plt.plot(s[1],s[3], marker = "o",color="r")
# for s in CE_sys:
#     plt.plot(s[1],s[3], marker = "o",color="y")
# plt.xlabel("Mass 1")
# plt.ylabel("Mass 2")
# # plt.legend(["Negative", "Positive", "CE"])
# plt.title("Before Mass Transfer/Common Envelope Formaiton")

# plt.subplot(2,1,2)
# for s in CEN_sys:
#     plt.plot(s[2],s[4], marker = "o", color = "b")
# for s in CEP_sys:
#     plt.plot(s[2],s[4], marker = "o",color="r")
# for s in CE_sys:
#     plt.plot(s[2],s[4], marker = "o",color="y")
# plt.xlabel("Mass 1")
# plt.ylabel("Mass 2")
# plt.legend(["Negative", "Positive", "CE"])
# plt.title("After Mass Transfer/Common Envelope Formaiton")
# plt.show()

# print(CEN_sys[0])
# print(CEP_sys[0])
# for CESeed in :
#     for seedIndex in range(len(seedsSP)):
#         spSeed = seedsSP[seedIndex]
#         if spSeed == dcSeed:
#             m1 = m1Zams[seedIndex]
#             m2 = m2Zams[seedIndex]
#             mTot = m1 + m2
#             totalMasses.append(mTot)


print(SPs.keys())
print(CEs.keys())

in_mass1 = SPs['Mass@ZAMS(1)'][()]
in_mass2 = SPs['Mass@ZAMS(2)'][()]
in_semax = SPs['SemiMajorAxis@ZAMS'][()]
#ce_alpha = SPs['CE_Alpha'][()]

seedsSP = SPs['SEED'][()]
seedsCE = CEs['SEED'][()]
# Calculate mask for which elements of seedsSP are found in seedsDC
# - see numpy.in1d documentation for details
mask = np.in1d(seedsSP, seedsCE)
# print(mask)
# print(seedsCE)
# print(seedsSP[mask])
print("The occurence rate of CEs is {}/{}".format(sum(mask), len(mask)))
print(len(seedsCE))
systems = np.arange(0, len(seeds),1)
plt.suptitle("Uniform IMF")
plt.subplot(4,1,1)
plt.hist([in_mass1, in_mass2, in_mass1+in_mass2])
plt.xlabel("Mass [Msun]")
plt.ylabel("# of Systems")
plt.yscale("log")
plt.legend(["Mp", "Ms", "Mtot"])

plt.subplot(4,1,2)
plt.hist(in_semax)
plt.xlabel("Semimajor axis [AU]")
plt.yscale("log")
plt.ylabel("# of Systems")

plt.subplot(4,1,3)
plt.hist([in_mass1[mask], CEs["Mass(1)<CE"],CEs["Mass(1)>CE"], in_mass2[mask],CEs["Mass(2)<CE"],CEs["Mass(2)>CE"] ])
plt.xlabel("Mass [Msun]")
plt.ylabel("CE positive")
plt.yscale("log")
plt.legend(["Mp_i", "Mp<CE", "Mp>CE","Ms_i", "Ms<CE", "Ms>CE"])

plt.subplot(4,1,4)
plt.hist([in_semax[mask], CEs["SemiMajorAxis<CE"], CEs["SemiMajorAxis>CE"]])
plt.xlabel("Semimajor axis [AU]")
plt.yscale("log")
plt.ylabel("CE Positive")
plt.legend(["a_i", "a<CE", "a>CE"])

plt.show()