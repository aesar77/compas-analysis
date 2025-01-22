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

models = ["Kroupa", "Salpeter", "Uniform", "Powerlaw"]
# plt.figure(figsize=(12, 12))

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with
for model in models:
    matplotlib.rcParams['figure.figsize'] = (15,10)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 14
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = '/data/a.saricaoglu/Runs/' + model + '/'
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

    mP_innms = []
    mP_ins = []
    mP_preCEs = []
    mP_pstCEs = []

    mS_innms = []
    mS_ins = []
    mS_preCEs = []
    mS_pstCEs = []

    semax_innms = []
    semax_ins = []
    semax_preCEs = []
    semax_pstCEs = []

    maskSP = np.array([], dtype='bool')
    # maskCEpre = np.array(0, dtype='bool')
    # maskCEpst = np.array(0, dtype='bool') #mask_unique
    maskCEMult = np.array([], dtype='bool')
    maskpreCEBH = np.array([], dtype='bool')
    maskpstCEBH = np.array([], dtype='bool')
    # maskpstCEBH =np.array(0, dtype='bool')
    maskpreCEBBH =np.array([], dtype='bool')
    maskpstCEBBH =np.array([], dtype='bool')

    maskSPBH = np.array([], dtype='bool')
    maskSPBBH =np.array([], dtype='bool')
    maskSPDCO = np.array([], dtype='bool')
    maskSPUnbound =  np.array([], dtype='bool')
    maskSPMerg = np.array([], dtype='bool')

    maskSPBH_CEpos = np.array([], dtype='bool')
    maskSPBBH_CEpos =np.array([], dtype='bool')
    maskSPDCO_CEpos = np.array([], dtype='bool')
    maskSPUnbound_CEpos =  np.array([], dtype='bool')
    maskSPMerg_CEpos = np.array([], dtype='bool')
    maskSPBoundBinary = np.array([], dtype='bool')
    maskSPBoundBinary_CEpos = np.array([], dtype='bool')

    # orbital_period = []
    # orbital_period_masked = []

    f = open("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Binary_BH/" + model + "_outputs.txt", "a")

    f.writelines(["\n","\n Run :", current_time])

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

        # print("SP KEYS: ", SP.keys())
        # print("CE KEYS: ", CE.keys())
        # print("DC KEYS: ",DC.keys())

    #IMPORTANT: a seed may appear more than once in CE[seeds], if they have more than one CE events. However, CE[CE_Event_Counter] output will have the same lenght with CE[seeds], and when
        #duplicates of the same seed appear in it as following, [..., 'X = 1', 'X = 2',.. ] so when we mask events = 1 we actually select all unique seeds, and for duplicates we consider the
        #state of the seed which it only had one CE event yet.
        seedsSP = np.unique(SP['SEED'][()])
        seedsCE = CE['SEED'][()]
        seedsDC = DC['SEED'][()]
        events   =  CE['CE_Event_Counter'][()]
        status = SP['Evolution_Status'][()]

        stellarTypeCE1   =  CE['Stellar_Type(1)'][()]
        stellarTypeCE2   =  CE['Stellar_Type(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]
        stellarTypeDC1   =  DC['Stellar_Type(1)'][()]
        stellarTypeDC2   =  DC['Stellar_Type(2)'][()]  

        maskCE  =  (events == 1)
        maskCEmult =  (events != 1) #masks all seeds having 1 events

        # print(len(CE['SEED'][()]), len(seedsSP), len(seedsDC))
        # print("len maskCE: ", len(maskCE), " sum maskCE: ", np.sum(maskCE))
        # print("len maskCEpst: ", len(maskCEpst), " sum maskCEpst: ", np.sum(maskCEpst))
        # print(maskCE)

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
        # print(mask_unique)
        # print(type(mask_unique), type((maskCE)) )
        # print("# of systems: ", len(seedsSP), "\n # of systems undergo CE: ", np.sum(maskCE), " \n $%$ of systems undergo CE: ",(np.sum(maskCE)*100)/len(seedsSP), 
        #       " \n # of systems undergo CE more than once: ",(np.sum(maskCEpst)), " \n $%$ of systems undergo CE: ",(np.sum(maskCEpst)*100/len(seedsSP)))
        
        #Selects systems which one of the pairs in the binary is a BH before CE formation
        maskCEBH1 = ((stellarTypeCE1[maskCE]  == 14) ^ (stellarTypeCE2[maskCE]  == 14))
        # #Selects systems which one of the pairs in the binary is a BH after CE formation
        maskCEBH2 = ((stellarTypeCE1[mask_unique]  == 14) ^ (stellarTypeCE2[mask_unique]  == 14)) 
        #Selects systems which both of the pairs in the binary are BH before CE formation
        maskCEBBH1 = ((stellarTypeCE1[maskCE]  == 14) & (stellarTypeCE2[maskCE]  == 14))
        #Selects systems which both of the pairs in the binary are BH after CE formation
        maskCEBBH2 = ((stellarTypeCE1[mask_unique]  == 14) & (stellarTypeCE2[mask_unique]  == 14))



        mask = np.in1d(seedsSP, seedsCE)
        print("pre CE mask: ", np.sum(mask), "post CE mask: ", np.sum(mask_unique))
        # mask_spdc = np.in1d(seedsSP, seedsDC)
        # mask_cedc = np.in1d(seedsCE, seedsDC)

        maskSPBH1 = ((stellarTypeSP1  == 14) ^ (stellarTypeSP2  == 14)) 
        maskSPBH2 = ((stellarTypeSP1  == 14) & (stellarTypeSP2  == 14))
        maskSPBH1_CEpos = ((stellarTypeSP1[mask]  == 14) ^ (stellarTypeSP2[mask]  == 14)) 
        maskSPBH2_CEpos= ((stellarTypeSP1[mask]  == 14) & (stellarTypeSP2[mask]  == 14)) 

        maskSPunb = (status == 14) & maskSPBH2
        maskSPdco = (status == 11) & maskSPBH2
        maskSPmrgr = (status == 9)  &  maskSPBH2
        maskSPunb_CEpos = (status[mask] == 14) & maskSPBH2_CEpos 
        maskSPdco_CEpos = (status[mask] == 11) & maskSPBH2_CEpos 
        maskSPmrgr_CEpos = (status[mask] == 9) & maskSPBH2_CEpos 
        #Note: if you were used status[mask] then maskSPunb_CEpos would have the same len with CE positive systems, and on the plots you would use mP_ins to mask.
        maskSPBndBnry = (maskSPBH2 == 1) & (maskSPdco == 1) & (maskSPunb == 0)
        maskSPBndBnry_CEpos = (maskSPBH2_CEpos == 1) & (maskSPdco_CEpos == 1) & (maskSPunb_CEpos == 0)
        print("pre CE maskpreCEBH: ", np.sum(maskCEBH1), "post CE maskpstCEBH: ", np.sum(maskCEBH2))


        # maskDCBH1 = ((stellarTypeDC1  == 14) | (stellarTypeDC2  == 14)) 
        # maskDCBH2 = ((stellarTypeDC1  == 14) & (stellarTypeDC2  == 14))
        # maskDCBH1_CEpos = ((stellarTypeDC1[mask]  == 14) | (stellarTypeDC2[mask]  == 14)) 
        # maskDCBH2_CEpos= ((stellarTypeDC1[mask]  == 14) & (stellarTypeDC2[mask]  == 14)) 

        # maskDCunb = (status == 14)
        # maskDCdco = (status == 11)
        # maskDCmrgr = (status == 9) 
        # maskDCunb_CEpos = (status[mask] == 14) 
        # maskDCdco_CEpos = (status[mask] == 11) 
        # maskDCmrgr_CEpos = (status[mask] == 9) 


        mP_innm = SP['Mass@ZAMS(1)'][()]
        mS_innm = SP['Mass@ZAMS(2)'][()]
        mP_in = SP['Mass@ZAMS(1)'][()][mask]
        mS_in = SP['Mass@ZAMS(2)'][()][mask]

        mP_preCE = CE["Mass(1)<CE"][()][maskCE]
        mS_preCE = CE["Mass(2)<CE"][()][maskCE]

        mP_pstCE = CE["Mass(1)>CE"][()][mask_unique]
        mS_pstCE = CE["Mass(2)>CE"][()][mask_unique]


        # mP_in = SP['Mass@ZAMS(1)'][()]
        # mS_in = SP['Mass@ZAMS(2)'][()]

        # mP_preCE = CE["Mass(1)<CE"][()][maskCE]
        # mS_preCE = CE["Mass(2)<CE"][()][maskCE]

        # mP_pstCE = CE["Mass(1)>CE"][()][mask_unique]
        # mS_pstCE = CE["Mass(2)>CE"][()][mask_unique]

        # mtP_in = SP['Mass@ZAMS(1)'][()][maskCEBH1]
        # mtS_in = SP['Mass@ZAMS(2)'][()][maskCEBH1]

        # mtP_preCE = CE["Mass(1)<CE"][()][maskCEBH1]
        # mtS_preCE = CE["Mass(2)<CE"][()][maskCEBH1]

        # mtP_pstCE = CE["Mass(1)>CE"][()][maskBH2]
        # mtS_pstCE = CE["Mass(2)>CE"][()][maskBH2]


        # orbital_period = SP['Orbital_Period'][()]
        semax_innm = SP['SemiMajorAxis@ZAMS'][()]
        semax_in = SP['SemiMajorAxis@ZAMS'][()][mask]
        semax_preCE = CE["SemiMajorAxis<CE"][()][maskCE]
        semax_pstCE = CE["SemiMajorAxis>CE"][()][mask_unique]

        mP_innms.extend(mP_innm)
        mP_ins.extend(mP_in)
        mP_preCEs.extend(mP_preCE)
        mP_pstCEs.extend(mP_pstCE)
        mS_innms.extend(mS_innm)
        mS_ins.extend(mS_in)
        mS_preCEs.extend(mS_preCE)
        mS_pstCEs.extend(mS_pstCE)

        # orbital_periods.extend(orbital_period)
        # orbital_periods_masked.extend(orbital_period[mask])
        semax_innms.extend(semax_innm)
        semax_ins.extend(semax_in)
        semax_preCEs.extend(semax_preCE)
        semax_pstCEs.extend(semax_pstCE)

        # np.concatenate([maskSP, mask)
        # np.concatenate([maskCEpre, maskCE)
        # np.concatenate([maskCEpst, mask_unique)
        maskpreCEBH = np.concatenate([maskpreCEBH, maskCEBH1])
        maskpstCEBH = np.concatenate([maskpstCEBH, maskCEBH2])
        maskpreCEBBH = np.concatenate([maskpreCEBBH, maskCEBBH1])
        maskpstCEBBH = np.concatenate([maskpstCEBBH, maskCEBBH2])
        maskSPBH = np.concatenate([maskSPBH, maskSPBH1])
        maskSPBBH = np.concatenate([maskSPBBH, maskSPBH2])
        maskSPUnbound = np.concatenate([maskSPUnbound, maskSPunb])
        maskSPDCO = np.concatenate([maskSPDCO, maskSPdco])
        maskSPMerg = np.concatenate([maskSPMerg, maskSPmrgr])
        maskSPBH_CEpos = np.concatenate([maskSPBH_CEpos, maskSPBH1_CEpos])
        maskSPBBH_CEpos = np.concatenate([maskSPBBH_CEpos, maskSPBH2_CEpos])
        maskSPUnbound_CEpos = np.concatenate([maskSPUnbound_CEpos, maskSPunb_CEpos])
        maskSPDCO_CEpos = np.concatenate([maskSPDCO_CEpos, maskSPdco_CEpos])
        maskSPMerg_CEpos = np.concatenate([maskSPMerg_CEpos, maskSPmrgr_CEpos])
        maskCEMult = np.concatenate([maskCEMult, maskCEmult])
        maskSPBoundBinary = np.concatenate([maskSPBoundBinary, maskSPBndBnry])
        maskSPBoundBinary_CEpos = np.concatenate([maskSPBoundBinary_CEpos, maskSPBndBnry_CEpos])

        Data.close()

    print(len(mP_ins), len(mP_preCEs), len(mS_pstCEs))
    print("BH percentage in CE positive systems: ", np.sum(maskpstCEBH)*100/len(maskpstCEBH), "\n BBH percentage in CE positive systems: ", np.sum(maskpstCEBBH)*100/len(maskpstCEBBH), 
        "\n BH percentage in all systems ", np.sum(maskSPBH)*100/len(maskSPBH), "\n BBH percentage in all systems: ", np.sum(maskSPBBH)*100/len(maskSPBH))

    L = [ "\n # of sys.: ",str(len(mP_innms)),  "\n # of sys. undergo CE: ", str(len(mP_ins)), "\n # of sys. undergo CE more than once: ",str(np.sum(maskCEMult)),  
        "\n $%$ of sys. undergo CE: ",str((len(mP_ins)*100)/len(mP_innms)),  " \n $%$ of sys. undergo CE more than once: ",str((np.sum(maskCEMult)*100/len(maskCEMult))), 
        "\n # of BH in all Systems: ", str(np.sum(maskSPBH)),  "\n # of (single) BH in CE pos sys.:  ", str(np.sum(maskSPBH_CEpos)),  
        "\n $%$ of BH in all sys.: ", str(np.sum(maskSPBH)*100/len(maskSPBH)),  "\n $%$ of (single) BH in CE pos sys.:", str(np.sum(maskSPBH_CEpos)*100/len(maskSPBH1_CEpos)),  
        "\n # of BBH in all Systems: ", str(np.sum(maskSPBBH)),  "\n # of BBH CE pos sys.:  ", str(np.sum(maskSPBBH_CEpos)),  
        "\n $%$ of BBH in all sys.:", str(np.sum(maskSPBBH)*100/len(maskSPBBH)),  "\n $%$ of BBH in CE pos sys.:", str(np.sum(maskSPBBH_CEpos)*100/len(maskSPBBH_CEpos)),  
        "\n # of unbound (2 BH) binaries in all sys.: ", str(np.sum(maskSPUnbound)), "\n # of unbound (2 BH) binaries in CE pos sys.: ", str(np.sum(maskSPUnbound_CEpos)), 
        "\n $%$ of unbound  (2 BH) binaries in all sys.: ", str(np.sum(maskSPUnbound)*100/len(maskSPUnbound)), "\n $%$ of unbound (2 BH) binaries in CE pos sys.: ", str(np.sum(maskSPUnbound_CEpos)*100/len(maskSPUnbound_CEpos)),  
        "\n # of double compact obj. (2 BH) in all sys.: ", str(np.sum(maskSPDCO)),  "\n # of double compact obj. (2 BH) in CE pos sys.: ", str(np.sum(maskSPDCO_CEpos)), 
        "\n $%$ of double compact obj. (2 BH) in all sys.: ", str(np.sum(maskSPDCO)*100/len(maskSPDCO)),  "\n $%$ of double compact obj. (2 BH) in CE pos sys.: ", str(np.sum(maskSPDCO_CEpos)*100/len(maskSPDCO_CEpos)), 
        "\n # of mergers (2 BH) in all sys.: ", str(np.sum(maskSPMerg)),  "\n # of (BBH) mergers in CE pos sys.: ", str(np.sum(maskSPMerg_CEpos)),  
        "\n $%$ of mergers (2 BH) in all sys.: ", str(np.sum(maskSPMerg)*100/len(maskSPMerg)), "\n $%$ of (2 BH) mergers in CE pos sys.: ", str(np.sum(maskSPMerg_CEpos)*100/len(maskSPMerg_CEpos)), 
        "\n # of bound binaries (2 BH) in all sys.: ", str(np.sum(maskSPBoundBinary)), "\n # of bound binaries (2 BH) in CE pos sys.: ", str(np.sum(maskSPBoundBinary_CEpos)),
        "\n $%$ of bound binaries (2 BH) in all sys.: ", str(np.sum(maskSPBoundBinary)*100/len(maskSPBoundBinary)), "\n $%$ of bound binaries (2 BH) in CE pos sys.: ", str(np.sum(maskSPBoundBinary_CEpos*100/len(maskSPBoundBinary_CEpos)))
        ]

    f.writelines(L)
    f.close()
    print("maskSPBH: ", maskSPBH, " maskSPBBH: ", maskSPBBH)

    #YOU INVERT THESE BECAUSE IN PLOTTING BELOW, ma.mask requires so
    maskSP = np.invert(maskSP)
    maskpstCEBH = np.invert(maskpstCEBH)
    maskpreCEBH = np.invert(maskpreCEBH)
    maskpreCEBBH = np.invert(maskpreCEBBH)
    maskpstCEBBH = np.invert(maskpstCEBBH)

    maskSPBH = np.invert(maskSPBH)
    maskSPBBH = np.invert(maskSPBBH)
    maskSPDCO = np.invert(maskSPDCO)
    maskSPUnbound = np.invert(maskSPUnbound)
    maskSPMerg = np.invert(maskSPMerg)
    maskSPBoundBinary = np.invert(maskSPBoundBinary)

    maskSPBH_CEpos = np.invert(maskSPBH_CEpos)
    maskSPBBH_CEpos = np.invert(maskSPBBH_CEpos)
    maskSPDCO_CEpos = np.invert(maskSPDCO_CEpos)
    maskSPUnbound_CEpos = np.invert(maskSPUnbound_CEpos)
    maskSPMerg_CEpos = np.invert(maskSPMerg_CEpos)
    maskSPBoundBinary_CEpos = np.invert(maskSPBoundBinary_CEpos)

    print("mP_ins: ", len(mP_ins), " mP_ins BH: ", len(ma.masked_array(mP_ins, mask=maskSPBH_CEpos)), " mP_ins BBH :", len(ma.masked_array(mP_ins, mask=maskSPBBH_CEpos)))
    print("mP_ins BH: ", np.sum(maskSPBH_CEpos), " mP_ins BBH :", np.sum(maskSPBBH_CEpos))

    # CEpos = []
    # CEneg = []
    # seeds, events = getEventHistory(Data)
    # print(events)

    # for ii, seed in enumerate([seeds):
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
    # # plt.close()

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



  
    # plt.subplot(3,3,1)
    # plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + ", CE positive: " + str(len(mP_ins)))
    # plt.hist([mP_innms, mP_ins], bins=40, log=True)
    # plt.xlabel("Mass [$M_\odot$]")
    # plt.ylabel("# of Systems")
    # plt.yscale("log")
    # #plt.xscale("log")
    # plt.legend(["All $M_p$ initial","$M_p$ initial undergoing CE"])

    # plt.subplot(3,3,2)
    # plt.suptitle(model+ " IMF")
    # plt.hist([mS_innms, mS_ins], bins=40, log=True)
    # plt.xlabel("Mass [$M_\odot$]")
    # plt.ylabel("# of Systems")
    # plt.yscale("log")
    # #plt.xscale("log")
    # plt.legend(["All $M_s$ initial","$M_s$ initial undergoing CE"])

    # plt.subplot(3,3,3)
    # plt.hist([semax_innms, semax_ins],bins=40, log=True)
    # plt.xlabel("Semimajor axis [AU]")
    # plt.yscale("log") 
    # #plt.xscale("log")
    # plt.ylabel("# of Systems")
    # plt.legend(["All $a$ initial","$a$ initial undergoing CE"])

    # plt.subplot(3,3,4)
    # plt.hist([mP_ins, mS_ins], bins=40, log=True)
    # plt.xlabel("Mass [$M_\odot$]")
    # plt.ylabel("CE positive Systems")
    # #plt.xscale("log")
    # plt.yscale("log")
    # plt.legend(["$M_p$ initial","$M_s$ initial"])

    # plt.subplot(3,3,5)
    # plt.hist([ mP_preCEs,mS_preCEs], bins=40, log=True)
    # plt.xlabel("Mass [$M_\odot$]")
    # plt.ylabel("CE positive Systems")
    # #plt.xscale("log")
    # plt.yscale("log")
    # plt.legend(["$M_p$ before common envelope", "$M_s$ before common envelope"])

    # plt.subplot(3,3,6)
    # plt.hist([mP_pstCEs,mS_pstCEs],  bins=40, log=True)
    # plt.xlabel("Mass [$M_\odot$]")
    # plt.ylabel("CE positive Systems")
    # #plt.xscale("log")
    # plt.yscale("log")
    # plt.legend(["$M_p$ after common envelope", "$M_s$ after common envelope"])

    # plt.subplot(3,3,7)
    # plt.hist([semax_ins, semax_preCEs],bins=40, log=True)
    # plt.xlabel("Semimajor axis [AU]")
    # plt.yscale("log")
    # #plt.xscale("log")
    # plt.ylabel("CE positive Systems")
    # plt.legend(["$a$ initial", "$a$ before common envelope"])

    # plt.subplot(3,3,8)
    # plt.hist([semax_preCEs, semax_pstCEs],bins=40, log=True)
    # plt.xlabel("Semimajor axis [AU]")
    # plt.yscale("log")
    # #plt.xscale("log")
    # plt.ylabel("CE positive Systems")
    # plt.legend(["$a$ before common envelope", "$a$ after common envelope"])

    # plt.subplot(3,3,9)
    # plt.hist([semax_ins, semax_preCEs],bins=40, log=True)
    # plt.xlabel("Semimajor axis [AU]")
    # plt.yscale("log")
    # #plt.xscale("log")
    # plt.ylabel("CE positive Systems")
    # plt.legend(["$a$ initial", "$a$ before common envelope"])



    # plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.3)
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + model + "_hist_" + current_time + ".png")
    # plt.close()

    x_scale = np.logspace(np.log10(1.0), np.log10(max(mP_ins)), 50)
    y_scale = np.logspace(np.log10(1.0), np.log10(max(mS_ins)), 50)

    plt.subplot(3,3,1)
    plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + "\n Evolution Status vs Masses in CE Positive Systems")
    plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPBH_CEpos)-ma.masked_array(mP_ins, mask= maskSPBBH_CEpos)), 
               abs(ma.masked_array(mS_ins, mask=maskSPBH_CEpos)-ma.masked_array(mS_ins, mask=maskSPBBH_CEpos)),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial  [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "(BH-BBH) Difference")

    plt.subplot(3,3,2)
    plt.hist2d(abs(ma.masked_array(mP_preCEs,mask=maskpreCEBH)-ma.masked_array(mP_preCEs, mask=maskpreCEBBH)), 
               abs(ma.masked_array(mS_preCEs, mask=maskpreCEBH)-ma.masked_array(mS_preCEs, mask=maskpreCEBBH)),bins=(x_scale, y_scale),)
    plt.xlabel("$M_p$ before CE [$M_{sun}$]")
    plt.ylabel("$M_s$ before CE [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "(BH-BBH) Difference")

    plt.subplot(3,3,3)
    plt.hist2d(abs(ma.masked_array(mP_pstCEs, mask=maskpstCEBH)-ma.masked_array(mP_pstCEs, mask=maskpstCEBBH)),
                (ma.masked_array(mS_pstCEs, mask=maskpstCEBH)-ma.masked_array(mS_pstCEs, mask=maskpstCEBBH)),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ after CE [$M_{sun}$]")
    plt.ylabel("$M_s$ after CE [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "(BH-BBH) Difference")

    plt.subplot(3,3,4)
    plt.hist2d(ma.masked_array(mP_ins, mask= maskSPBH_CEpos), ma.masked_array(mS_ins, mask=maskSPBH_CEpos),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial CE [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "Black Hole occurance")

    plt.subplot(3,3,5)
    plt.hist2d(ma.masked_array(mP_ins, mask= maskSPBBH_CEpos), ma.masked_array(mS_ins, mask=maskSPBBH_CEpos),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial CE [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "Binary Black Hole occurance")

    plt.subplot(3,3,6)
    plt.hist2d(ma.masked_array(mP_ins, mask=maskSPUnbound_CEpos), ma.masked_array(mS_ins, mask=maskSPUnbound_CEpos),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "Unbound Binary (1 BH) occurance")

    plt.subplot(3,3,7)
    plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPDCO_CEpos)-ma.masked_array(mP_ins, mask= maskSPBBH_CEpos)),
                abs(ma.masked_array(mS_ins, mask=maskSPDCO_CEpos)-ma.masked_array(mS_ins, mask= maskSPBBH_CEpos)),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "(DCO-BBH) Difference")

    plt.subplot(3,3,8)
    plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPUnbound_CEpos)-ma.masked_array(mP_ins, mask=maskSPMerg_CEpos)), 
               abs(ma.masked_array(mS_ins, mask=maskSPUnbound_CEpos)- ma.masked_array(mS_ins, mask=maskSPMerg_CEpos)),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial[$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "(UBB-Merg) Difference")

    plt.subplot(3,3,9)
    plt.hist2d(ma.masked_array(mP_ins, mask=maskSPDCO_CEpos), ma.masked_array(mS_ins, mask=maskSPDCO_CEpos),bins=(x_scale, y_scale))
    plt.xlabel("$M_p$ initial [$M_{sun}$]")
    plt.ylabel("$M_s$ initial [$M_{sun}$]")
    plt.yscale("log")
    plt.xscale("log")
    plt.colorbar(label = "DCO (except BBH) occurance")


    plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.45)
    plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Binary_BH/" + model + "_2dhist_" + current_time + ".png")
    plt.close()

    plt.subplot(3,3,1)
    plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + ", CE positive: " + str(len(mP_ins)))
    plt.hist([mP_innms, mP_ins], bins=40, log=True)
    plt.xlabel("Mass [$M_\odot$]")
    plt.ylabel("# of Systems")
    plt.yscale("log")
    #plt.xscale("log")
    plt.legend(["All $M_p$ initial","$M_p$ initial undergoing CE"])

    plt.subplot(3,3,2)
    plt.suptitle(model+ " IMF")
    plt.hist([mS_innms, mS_ins], bins=40, log=True)
    plt.xlabel("Mass [$M_\odot$]")
    plt.ylabel("# of Systems")
    plt.yscale("log")
    #plt.xscale("log")
    plt.legend(["All $M_s$ initial","$M_s$ initial undergoing CE"])

    plt.subplot(3,3,3)
    plt.hist([semax_innms, semax_ins],bins=40, log=True)
    plt.xlabel("Semimajor axis [AU]")
    plt.yscale("log") 
    #plt.xscale("log")
    plt.ylabel("# of Systems")
    plt.legend(["All $a$ initial","$a$ initial undergoing CE"])

    plt.subplot(3,3,4)
    plt.hist([mP_ins, mS_ins], bins=40, log=True)
    plt.xlabel("Mass [$M_\odot$]")
    plt.ylabel("CE positive Systems")
    #plt.xscale("log")
    plt.yscale("log")
    plt.legend(["$M_p$ initial","$M_s$ initial"])

    plt.subplot(3,3,5)
    plt.hist([ mP_preCEs,mS_preCEs], bins=40, log=True)
    plt.xlabel("Mass [$M_\odot$]")
    plt.ylabel("CE positive Systems")
    #plt.xscale("log")
    plt.yscale("log")
    plt.legend(["$M_p$ before common envelope", "$M_s$ before common envelope"])

    plt.subplot(3,3,6)
    plt.hist([mP_pstCEs,mS_pstCEs],  bins=40, log=True)
    plt.xlabel("Mass [$M_\odot$]")
    plt.ylabel("CE positive Systems")
    #plt.xscale("log")
    plt.yscale("log")
    plt.legend(["$M_p$ after common envelope", "$M_s$ after common envelope"])

    plt.subplot(3,3,7)
    plt.hist([semax_ins, semax_preCEs],bins=40, log=True)
    plt.xlabel("Semimajor axis [AU]")
    plt.yscale("log")
    #plt.xscale("log")
    plt.ylabel("CE positive Systems")
    plt.legend(["$a$ initial", "$a$ before common envelope"])

    plt.subplot(3,3,8)
    plt.hist([semax_preCEs, semax_pstCEs],bins=40, log=True)
    plt.xlabel("Semimajor axis [AU]")
    plt.yscale("log")
    #plt.xscale("log")
    plt.ylabel("CE positive Systems")
    plt.legend(["$a$ before common envelope", "$a$ after common envelope"])

    plt.subplot(3,3,9)
    plt.hist([semax_ins, semax_preCEs],bins=40, log=True)
    plt.xlabel("Semimajor axis [AU]")
    plt.yscale("log")
    #plt.xscale("log")
    plt.ylabel("CE positive Systems")
    plt.legend(["$a$ initial", "$a$ before common envelope"])



    plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.3)
    plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Binary_BH/" + model + "_hist_" + current_time + ".png")
    plt.close()

    y_scale2 = np.logspace(np.log10(1.0), np.log10(max(semax_innms)), 50)



    plt.subplot(3,2,1)
    plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + "\n Evolution Status vs Semimajor Axis in CE Positive Systems")
    plt.hist([ma.masked_array(semax_ins, mask=maskSPBH_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" Single BH Binaries ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplot(3,2,2)
    plt.hist([ma.masked_array(semax_ins, mask=maskSPBBH_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" BH Binaries ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplot(3,2,3)
    plt.hist([ma.masked_array(semax_ins, mask=maskSPBoundBinary_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" Bound Single BH Binaries ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplot(3,2,4)
    plt.hist([ma.masked_array(semax_ins, mask=maskSPDCO_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" DCOs w/ Single BH ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplot(3,2,5)
    plt.hist([ma.masked_array(semax_ins, mask=maskSPUnbound_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" Unbound Sys. w/ Single BH ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplot(3,2,6)
    plt.hist([ma.masked_array(semax_ins, mask=maskSPMerg_CEpos)],bins=50)
    plt.xlabel("$a$ initial [AU]")
    plt.ylabel(" Mergers w/ Single BH ")
    plt.yscale("log")
    #plt.xscale("log")

    plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.45)
    plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Binary_BH/" + model + "_histsemax_" + current_time + ".png")
    plt.close()
    # plt.subplot(3,3,1)
    # plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + ", CE positive: " + str(len(mP_ins)))
    # plt.scatter(semax_innms, mP_innms, marker='.')
    # plt.scatter(semax_ins, mP_ins ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_p$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["All Systems","Systems undergoing CE"])

    # plt.subplot(3,3,2)
    # plt.scatter(semax_innms, mS_innms ,marker='.')
    # plt.scatter(semax_ins, mP_ins ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["All Systems","Systems undergoing CE"])

    # plt.subplot(3,3,3)
    # plt.scatter(semax_innms, np.add(mS_innms, mP_innms) ,marker='.')
    # plt.scatter(semax_ins, np.add(mS_ins, mP_ins) ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s + M_p$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["All Systems","Systems undergoing CE"])


    # plt.subplot(3,3,4)
    # plt.scatter(semax_preCEs, mP_preCEs ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_p$ before CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["Systems undergoing CE"])

    # plt.subplot(3,3,5)
    # plt.scatter(semax_preCEs, mS_preCEs ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s$ before CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["Systems undergoing CE"])

    # plt.subplot(3,3,6)
    # plt.scatter(semax_preCEs, np.add(mS_preCEs,mS_preCEs) ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s + M_p$ before CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["Systems undergoing CE"])

    # plt.subplot(3,3,7)
    # plt.scatter(semax_pstCEs, mP_pstCEs,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_p$ after CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")

    # plt.legend(["Systems undergoing CE"])


    # plt.subplot(3,3,8)
    # plt.scatter(semax_pstCEs, mS_pstCEs,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s$ after CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.legend(["Systems undergoing CE"])

    # plt.subplot(3,3,9)
    # plt.scatter(semax_pstCEs, np.add(mS_pstCEs, mP_pstCEs) ,marker='.')
    # plt.xlabel("Semimajor Axis [AU]")
    # plt.ylabel("$M_s + M_p$ after CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.ylim(0,50)
    # plt.xlim(0,50)
    # plt.legend(["Systems undergoing CE"])


    # plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.235)
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + model +"_scatter_" + current_time + ".png")
    # plt.close()


    # matplotlib.rcParams['figure.figsize'] = (6,5)
    # matplotlib.rcParams['font.size'] = 10
    # plt.title("Mass Distribution for Systems with Bound Binaries")
    # plt.scatter(ma.masked_array(mP_innms, mask=maskSPBoundBinary),ma.masked_array(mS_innms, mask=maskSPBoundBinary) ,marker='.')
    # plt.scatter(ma.masked_array(mP_ins, mask=maskSPBoundBinary_CEpos),ma.masked_array(mS_ins, mask=maskSPBoundBinary_CEpos) ,marker='.')
    # plt.xlabel("$M_p$ (BH component) [$M_{sun}$]")
    # plt.ylabel("$M_s$ (BH/CO/Star component) [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # # plt.axvline(0.1, label="0.1 [$M_{sun}$]")
    # # plt.axvline(1000, label="1000 [$M_{sun}$]")
    # plt.legend(["All systems", "Systems undergoing CE"])
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + model + "_boundbh_" + current_time + ".png")
    # plt.close()

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
    # plt.close()

