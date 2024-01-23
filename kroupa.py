import os, sys
import numpy as np               # for handling arrays
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib

c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)

# plt.figure(figsize=(12, 12))
matplotlib.rcParams['figure.figsize'] = (15,10)
matplotlib.rcParams['lines.markersize'] = 1
matplotlib.rcParams['font.size'] = 14

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with

pathToData = '/data/a.saricaoglu/Runs/Kroupa/'
runs= [x[0] for x in os.walk(pathToData) if "COMPAS" in x[0]]
data_outputs = []
for run in runs:
    out = [f for f in os.listdir(run) if ".h5" in f]
    try:
        data = h5.File(run + "/" + out[0])
        print("Reading file: ", run + "/" + out[0])
        data_outputs.append(data)
    except:
        continue

SPs = []
MTs = []
CEs = []
SNe = []
DCs = []

mP_ins = []
mP_ins_masked = []
mP_preCEs = []
mP_pstCEs = []
mS_ins = []
mS_ins_masked = []
mS_preCEs = []
mS_pstCEs = []

mtP_ins = []
mtP_preCEs = []
mtP_pstCEs = []
mtS_ins = []
mtS_preCEs = []
mtS_pstCEs = []

semax_ins = []
semax_ins_masked = []
semax_preCEs = []
semax_pstCEs = []
# orbital_period = []
# orbital_period_masked = []
 
for Data in data_outputs:

    SP = Data['BSE_System_Parameters']
    # MT = Data['BSE_RLOF']
    CE = Data['BSE_Common_Envelopes']
    # SN = Data['BSE_Supernovae']
    DC = Data['BSE_Double_Compact_Objects']

    SPs.extend(SP)
    #MTs.extend(MT)
    CEs.extend(CE)
    #SNe.extend(SN)
    DCs.extend(DC)

    print(SP.keys())
    print(DC.keys())

#IMPORTANT: a seed may appear more than once in CE[seeds], if they have more than one CE events. However, CE[CE_Event_Counter] output will have the same lenght with CE[seeds], and when
    #duplicates of the same seed appear in it as following, [..., 'X = 1', 'X = 2',.. ] so when we mask events = 1 we actually select all unique seeds, and for duplicates we consider the
    #state of the seed which it only had one CE event yet.
    seedsSP = np.unique(SP['SEED'][()])
    seedsCE = CE['SEED'][()]
    seedsDC = DC['SEED'][()]
    events   =  CE['CE_Event_Counter'][()]
    status =  CE['Evolution_Status'][()]
    stellarTypeCE1   =  CE['Stellar_Type(1)'][()]
    stellarTypeCE2   =  CE['Stellar_Type(2)'][()]
    maskCE  =  (events == 1)
    maskCE2 =  (events != 1) #masks all seeds having 1 events
    maskCE3 = (status == 14)
    maskCE4 = (status == 11)
    #seedsCE = seedsCE[maskCE]
    print(len(CE['SEED'][()]), len(seedsSP), len(seedsDC))
    print("len maskCE: ", len(maskCE), " sum maskCE: ", np.sum(maskCE))
    print("len maskCE2: ", len(maskCE2), " sum maskCE2: ", np.sum(maskCE2))
    print(maskCE)

    ##Constructs mask_unique in a way to select unique seeds from CE[seeds], if a seed has duplicates in CE[seeds] then it chooses the first one (same as maskCE)
    # pre = 0
    # mask_unique = []
    # for i in range(0, len(CE['SEED'][()])):
    #     if (CE['SEED'][()][i] == pre):
    #         mask_unique.append(0)
    #        
    #         #print(events[i-1], events[i])
    #     else:
    #         mask_unique.append(1)
    #         pre = CE['SEED'][()][i]

#IMPORTANT: mask_unique is constructed in a way to select unique seeds from CE[seeds], if a seed has duplicates in CE[seeds] then it chooses the latest one (the one with the most CE event count)
    pre = 0
    mask_unique = np.zeros(len(CE['SEED'][()]), dtype='bool')
    for i in range(0, len(CE['SEED'][()])):
        if (CE['SEED'][()][i] == pre):
            mask_unique[i-1] = False
            mask_unique[i] = True
            #
            #print(events[i-1], events[i])
        else:
            mask_unique[i] = True
            pre = CE['SEED'][()][i]
    print(mask_unique)
    print(type(mask_unique), type((maskCE)) )
    print("Number of systems: ", len(seedsSP), "\n Number of systems undergo CE: ", np.sum(maskCE), " \n Percentage of systems undergo CE: ",(np.sum(maskCE)*100)/len(seedsSP), 
          " \n Number of systems undergo CE more than once: ",(np.sum(maskCE2)), " \n Percentage of systems undergo CE: ",(np.sum(maskCE2)*100/len(seedsSP)))
    
    #Selects systems which one of the pairs in the binary is a BH before CE formation
    maskBH1 = ((stellarTypeCE1  == 14) | (stellarTypeCE2  == 14)) & (maskCE == 1)
    #Selects systems which one of the pairs in the binary is a BH after CE formation
    maskBH2 = ((stellarTypeCE1  == 14) | (stellarTypeCE2  == 14)) & (mask_unique == 1)
    #Selects systems which both of the pairs in the binary are BH after CE formation
    maskBH3 = ((stellarTypeCE1  == 14) & (stellarTypeCE2  == 14)) & (mask_unique == 1)


    print(np.sum(maskBH1),np.sum(maskBH2))

    mP_in = SP['Mass@ZAMS(1)'][()]
    mS_in = SP['Mass@ZAMS(2)'][()]

    mP_preCE = CE["Mass(1)<CE"][()][maskCE]
    mS_preCE = CE["Mass(2)<CE"][()][maskCE]

    mP_pstCE = CE["Mass(1)>CE"][()][mask_unique]
    mS_pstCE = CE["Mass(2)>CE"][()][mask_unique]

    mtP_in = SP['Mass@ZAMS(1)'][()][maskBH1]
    mtS_in = SP['Mass@ZAMS(2)'][()][maskBH1]

    mtP_preCE = CE["Mass(1)<CE"][()][maskBH1]
    mtS_preCE = CE["Mass(2)<CE"][()][maskBH1]

    mtP_pstCE = CE["Mass(1)>CE"][()][maskBH2]
    mtS_pstCE = CE["Mass(2)>CE"][()][maskBH2]

    mask = np.in1d(seedsSP, seedsCE[maskCE])

    # orbital_period = SP['Orbital_Period'][()]

    semax_in = SP['SemiMajorAxis@ZAMS'][()]
    semax_preCE = CE["SemiMajorAxis<CE"][()][maskCE]
    semax_pstCE = CE["SemiMajorAxis>CE"][()][mask_unique]

    mP_ins.extend(mP_in)
    mP_preCEs.extend(mP_preCE)
    mP_pstCEs.extend(mP_pstCE)
    mS_ins.extend(mS_in)
    mS_preCEs.extend(mS_preCE)
    mS_pstCEs.extend(mS_pstCE)

    mtP_ins.extend(mtP_in)
    mtP_preCEs.extend(mtP_preCE)
    mtP_pstCEs.extend(mtP_pstCE)
    mtS_ins.extend(mtS_in)
    mtS_preCEs.extend(mtS_preCE)
    mtS_pstCEs.extend(mtS_pstCE)
    # orbital_periods.extend(orbital_period)
    # orbital_periods_masked.extend(orbital_period[mask])
    semax_ins.extend(semax_in)
    semax_ins_masked.extend(semax_in[mask])
    semax_preCEs.extend(semax_preCE)
    semax_pstCEs.extend(semax_pstCE)

    Data.close()

print(len(mP_ins), len(mP_ins_masked), len(mP_preCEs), len(mS_pstCEs))



# CEpos = []
# CEneg = []

# seeds, events = getEventHistory(Data)
# print(events)

# for ii, seed in enumerate(seeds):
#     # print(seed, events[ii])
#     occurance = False
#     for event in events[ii]:

#         if event[0] == "MT" and event[-1] == True and occurance == False:
#             #print(event[0], event[-1])
#             CEpos.append(seed)
#             occurance = True
#         elif seed not in CEneg:
#             CEneg.append(seed)

# # print(len(CEpos))
# # print(len(CEneg))
# # print(len(seeds))

# CEP_sys = []
# CEN_sys = []
# CE_sys = []

# # for negseed in CEneg:
# #     for mtSeedIndex in range(len(MTs['SEED'][()])):
# #         mtSeed = MTs['SEED'][mtSeedIndex]
# #         if mtSeed == negseed:
# #             CEN_sys.append([negseed, MTs["Mass(1)<MT"][mtSeedIndex],MTs["Mass(1)>MT"][mtSeedIndex],MTs["Mass(2)<MT"][mtSeedIndex],MTs["Mass(2)>MT"][mtSeedIndex], MTs['Eccentricity<MT'][mtSeedIndex], MTs['Eccentricity>MT'][mtSeedIndex]])

# # for posseed in CEpos:
# #     for ceSeedIndex in range(len(CEs['SEED'][()])):
# #         ceSeed = CEs['SEED'][ceSeedIndex]
# #         if ceSeed == posseed:
# #             CE_sys.append([posseed, CEs["Mass(1)<CE"][ceSeedIndex],CEs["Mass(1)>CE"][ceSeedIndex],CEs["Mass(2)<CE"][ceSeedIndex],CEs["Mass(2)>CE"][ceSeedIndex], CEs['Eccentricity<CE'][ceSeedIndex], CEs['Eccentricity>CE'][ceSeedIndex]])
# #     for mtSeedIndex in range(len(MTs['SEED'][()])):
# #         mtSeed = MTs['SEED'][mtSeedIndex]
# #         if mtSeed == posseed:
# #             CEP_sys.append([posseed, MTs["Mass(1)<MT"][mtSeedIndex],MTs["Mass(1)>MT"][mtSeedIndex],MTs["Mass(2)<MT"][mtSeedIndex],MTs["Mass(2)>MT"][mtSeedIndex], MTs['Eccentricity<MT'][mtSeedIndex], MTs['Eccentricity>MT'][mtSeedIndex]])

# # plt.subplot(2,1,1)
# # for s in CEN_sys:
# #     plt.plot(s[1],s[3], marker = "o", color = "b")
# # for s in CEP_sys:
# #     plt.plot(s[1],s[3], marker = "o",color="r")
# # for s in CE_sys:
# #     plt.plot(s[1],s[3], marker = "o",color="y")
# # plt.xlabel("Mass 1")
# # plt.ylabel("Mass 2")
# # # plt.legend(["Negative", "Positive", "CE"])
# # plt.title("Before Mass Transfer/Common Envelope Formaiton")

# # plt.subplot(2,1,2)
# # for s in CEN_sys:
# #     plt.plot(s[2],s[4], marker = "o", color = "b")
# # for s in CEP_sys:
# #     plt.plot(s[2],s[4], marker = "o",color="r")
# # for s in CE_sys:
# #     plt.plot(s[2],s[4], marker = "o",color="y")
# # plt.xlabel("Mass 1")
# # plt.ylabel("Mass 2")
# # plt.legend(["Negative", "Positive", "CE"])
# # plt.title("After Mass Transfer/Common Envelope Formaiton")
# # plt.show()

# # print(CEN_sys[0])
# # print(CEP_sys[0])
# # for CESeed in :
# #     for seedIndex in range(len(seedsSP)):
# #         spSeed = seedsSP[seedIndex]
# #         if spSeed == dcSeed:
# #             m1 = m1Zams[seedIndex]
# #             m2 = m2Zams[seedIndex]
# #             mTot = m1 + m2
# #             totalMasses.append(mTot)


# print(SPs.keys())
# print(CEs.keys())

# seedsSP = SPs['SEED'][()]
# seedsCE = np.unique(CEs['SEED'][()])

# mP_in = SPs['Mass@ZAMS(1)'][()]
# mS_in = SPs['Mass@ZAMS(2)'][()]

# mP_preCE = CEs["Mass(1)<CE"][()]
# mS_preCE = CEs["Mass(2)<CE"][()]

# mP_pstCE = CEs["Mass(1)>CE"][()]
# mS_pstCE = CEs["Mass(2)>CE"][()]

# mask_mp = np.in1d(seedsSP, seedsCE)
# mask_ms = np.in1d(seedsSP, seedsCE)

# semax_in = SPs['SemiMajorAxis@ZAMS'][()]
# semax_preCE = CEs["SemiMajorAxis<CE"]
# semax_pstCE = CEs["SemiMajorAxis>CE"]

# mask_a = np.in1d(seedsSP, seedsCE)



# # Calculate mask for which elements of seedsSP are found in seedsDC
# # - see numpy.in1d documentation for details
# mask = np.in1d(seedsSP, seedsCE)
# maskedseeds = np.in1d(seedsCE, seedsSP[mask])
# # for ceSeedIndex in range(len(CEs['SEED'][()])):
# #     for ceeSeed in seedsSP[mask]:
# #         ceSeed = CEs['SEED'][ceSeedIndex]
# #         if ceeSeed == ceSeed:
# #             CE_sys.append([ceeSeed, CEs["Mass(1)<CE"][ceSeedIndex],CEs["Mass(1)>CE"][ceSeedIndex],CEs["Mass(2)<CE"][ceSeedIndex],CEs["Mass(2)>CE"][ceSeedIndex], CEs['SemiMajorAxis<CE'][ceSeedIndex], CEs['SemiMajorAxis>CE'][ceSeedIndex]])

# # print(mask)
# # print(seedsCE)
# # print(seedsSP[mask])
# print("The occurence rate of CEs is {}/{}".format(sum(mask), len(mask)))
# print(len(np.in1d(mP_in[mask_mp], mP_preCE)))
# print(len(np.unique(seedsCE)))
# systems = np.arange(0, len(seeds),1)


plt.subplot(3,3,1)
plt.suptitle("Kroupa IMF, Population size: " + str(len((mP_ins))) + ", CE positive: " + str(len(mP_ins_masked)))
plt.hist([mP_ins, mP_ins_masked], bins=40)
plt.xlabel("Mass [$M_\odot$]")
plt.ylabel("# of Systems")
plt.yscale("log")
#plt.xscale("log")
plt.legend(["All $M_p$ initial","$M_p$ initial undergoing CE"])

plt.subplot(3,3,2)
plt.suptitle("Kroupa IMF")
plt.hist([mS_ins, mS_ins_masked], bins=40)
plt.xlabel("Mass [$M_\odot$]")
plt.ylabel("# of Systems")
plt.yscale("log")
#plt.xscale("log")
plt.legend(["All $M_s$ initial","$M_s$ initial undergoing CE"])

plt.subplot(3,3,3)
plt.hist([semax_ins, semax_ins_masked], bins=40)
plt.xlabel("Semimajor axis [AU]")
plt.yscale("log") 
#plt.xscale("log")
plt.ylabel("# of Systems")
plt.legend(["All $a$ initial","$a$ initial undergoing CE"])

plt.subplot(3,3,4)
plt.hist([mP_ins_masked, mS_ins_masked],  bins=40 )
plt.xlabel("Mass [$M_\odot$]")
plt.ylabel("CE positive Systems")
#plt.xscale("log")
plt.yscale("log")
plt.legend(["$M_p$ initial","$M_s$ initial"])

plt.subplot(3,3,5)
plt.hist([ mP_preCEs,mS_preCEs],  bins=40 )
plt.xlabel("Mass [$M_\odot$]")
plt.ylabel("CE positive Systems")
#plt.xscale("log")
plt.yscale("log")
plt.legend(["$M_p$ before common envelope", "$M_s$ before common envelope"])

plt.subplot(3,3,6)
plt.hist([mP_pstCEs,mS_pstCEs],  bins=40 )
plt.xlabel("Mass [$M_\odot$]")
plt.ylabel("CE positive Systems")
#plt.xscale("log")
plt.yscale("log")
plt.legend(["$M_p$ after common envelope", "$M_s$ after common envelope"])

plt.subplot(3,3,7)
plt.hist([semax_preCEs, semax_pstCEs], bins=40)
plt.xlabel("Semimajor axis [AU]")
plt.yscale("log")
#plt.xscale("log")
plt.ylabel("CE positive Systems")
plt.legend(["$a$ initial", "$a$ before common envelope", "$a$ after common envelope"])

plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.235)
plt.savefig("/data/a.saricaoglu/Plots/kroupa_hist_" + current_time + ".png")
plt.show()

plt.subplot(3,3,1)
plt.hist2d(mtP_in, mtS_in)
plt.xlabel("$M_p$ before CE [$M_{sun}$]")
plt.ylabel("$M_s$ before CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")

plt.subplot(3,3,2)
plt.hist2d(mtP_preCE, mtS_preCE)
plt.xlabel("$M_p$ before CE [$M_{sun}$]")
plt.ylabel("$M_s$ before CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")

plt.subplot(3,3,3)
plt.hist2d(mtP_pstCE, mtS_pstCE)
plt.xlabel("$M_p$ after CE [$M_{sun}$]")
plt.ylabel("$M_s$ after CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")

plt.subplot(3,3,4)
plt.hist2d(mtP_in, mtS_in)
plt.xlabel("$M_p$ before CE [$M_{sun}$]")
plt.ylabel("$M_s$ before CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")

plt.subplot(3,3,5)
plt.hist2d(mtP_preCE, mtS_preCE)
plt.xlabel("$M_p$ before CE [$M_{sun}$]")
plt.ylabel("$M_s$ before CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")

plt.subplot(3,3,6)
plt.hist2d(mtP_pstCE, mtS_pstCE)
plt.xlabel("$M_p$ after CE [$M_{sun}$]")
plt.ylabel("$M_s$ after CE [$M_{sun}$]")
#plt.yscale("log")
#plt.xscale("log")
plt.colorbar(label = "Black Hole occurance")



plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.235)
plt.savefig("/data/a.saricaoglu/Plots/kroupa_2dhist" + current_time + ".png")
plt.show()

plt.subplot(3,3,1)
plt.suptitle("Kroupa IMF, Population size: " + str(len((mP_ins))) + ", CE positive: " + str(len(mP_ins_masked)))
plt.scatter(semax_ins, mP_ins, marker='.')
plt.scatter(semax_ins_masked, mP_ins_masked ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_p$ initial [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["All Systems","Systems undergoing CE"])

plt.subplot(3,3,2)
plt.scatter(semax_ins, mS_ins ,marker='.')
plt.scatter(semax_ins_masked, mP_ins_masked ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s$ initial [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["All Systems","Systems undergoing CE"])

plt.subplot(3,3,3)
plt.scatter(semax_ins, np.add(mS_ins, mP_ins) ,marker='.')
plt.scatter(semax_ins_masked, mP_ins_masked ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s + M_p$ initial [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["All Systems","Systems undergoing CE"])


plt.subplot(3,3,4)
plt.scatter(semax_preCEs, mP_preCEs ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_p$ before CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["Systems undergoing CE"])

plt.subplot(3,3,5)
plt.scatter(semax_preCEs, mS_preCEs ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s$ before CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["Systems undergoing CE"])

plt.subplot(3,3,6)
plt.scatter(semax_preCEs, np.add(mS_preCEs,mS_preCEs) ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s + M_p$ before CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["Systems undergoing CE"])

plt.subplot(3,3,7)
plt.scatter(semax_pstCEs, mP_pstCEs,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_p$ after CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")

plt.legend(["Systems undergoing CE"])


plt.subplot(3,3,8)
plt.scatter(semax_pstCEs, mS_pstCEs,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s$ after CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["Systems undergoing CE"])

plt.subplot(3,3,9)
plt.scatter(semax_pstCEs, np.add(mS_pstCEs, mP_pstCEs) ,marker='.')
plt.xlabel("Semimajor Axis [AU]")
plt.ylabel("$M_s + M_p$ after CE [$M_{sun}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend(["Systems undergoing CE"])


plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.235)
plt.savefig("/data/a.saricaoglu/Plots/kroupa_scatter" + current_time + ".png")
plt.show()

# plt.subplot(1,1)
# plt.scatter(orbital_period, mP_ins, marker='.')
# plt.scatter(orbital_period_masked, mP_ins_masked, marker='.')
# plt.xlabel("Orbital Period [days]")
# plt.ylabel("$M_p$ initial [$M_{sun}$]")
# plt.yscale("log")
# #plt.xscale("log")
# plt.legend(["All systems", "Systems undergoing CE"])

# plt.subplot(1,2)
# plt.scatter(orbital_period, mS_ins, marker='.')
# plt.scatter(orbital_period_masked, mS_ins_masked, marker='.')
# plt.xlabel("Orbital Period [days]")
# plt.ylabel("$M_s$ initial [$M_{sun}$]")
# plt.yscale("log")
# #plt.xscale("log")
# plt.legend(["All systems", "Systems undergoing CE"])

# plt.subplots_adjust(top=0.93, bottom=0.07, left=0.045, right=0.955, hspace=0.32, wspace=0.235)
# plt.savefig("/data/a.saricaoglu/Plots/kroupa_T" + current_time + ".png")
# plt.show()

