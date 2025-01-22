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
print("current time : ", current_time)

models = ["Kroupa"]
mode = "SingleBH"
# plt.figure(figsize=(12, 12))

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with
for model in models:

    f = open("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode + "_" + model + "_outputs.txt", "a")

    f.writelines(["\n","\n Run : ", current_time])

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
    mP_DCs = []
    mP_DCs_CEs = []

    mS_innms = []
    mS_ins = []
    mS_preCEs = []
    mS_pstCEs = []
    mS_DCs = []
    mS_DCs_CEs =[]

    semax_innms = []
    semax_ins = []
    semax_preCEs = []
    semax_pstCEs = []
    semax_DCs = []
    semax_DCs_CEs = []

    final_status = np.array([], dtype='bool')

    maskSP = np.array([], dtype='bool')
    maskSPCE = np.array([], dtype='bool')

    # maskCEpre = np.array(0, dtype='bool')
    # maskCEpst = np.array(0, dtype='bool') #mask_unique
    maskCEMult = np.array([], dtype='bool')
    maskpreCEBH = np.array([], dtype='bool')
    maskpstCEBH = np.array([], dtype='bool')
    # maskpstCEBH =np.array(0, dtype='bool')
    maskpreCEBBH =np.array([], dtype='bool')
    maskpstCEBBH =np.array([], dtype='bool')

    maskSPBH = np.array([], dtype='bool')
    maskSPBHpr = np.array([], dtype='bool')
    maskSPBHsc = np.array([], dtype='bool')

    maskDCBH = np.array([], dtype='bool')
    maskDCBHpr = np.array([], dtype='bool')
    maskDCBHsc = np.array([], dtype='bool')

    maskSPBBH =np.array([], dtype='bool')
    maskDCBBH =np.array([], dtype='bool')


    maskSPDCO = np.array([], dtype='bool')
    maskSPUnbound =  np.array([], dtype='bool')
    maskSPMerg = np.array([], dtype='bool')

    maskSPBH_CEpos = np.array([], dtype='bool')
    maskSPBHpr_CEpos = np.array([], dtype='bool')
    maskSPBHsc_CEpos = np.array([], dtype='bool')

    maskDCBH_CEpos = np.array([], dtype='bool')
    maskDCBHpr_CEpos = np.array([], dtype='bool')
    maskDCBHsc_CEpos = np.array([], dtype='bool')

    maskSPBBH_CEpos =np.array([], dtype='bool')
    maskDCBBH_CEpos =np.array([], dtype='bool')


    maskSPDCO_CEpos = np.array([], dtype='bool')
    maskSPUnbound_CEpos =  np.array([], dtype='bool')
    maskSPMerg_CEpos = np.array([], dtype='bool')
    maskSPBoundBinary = np.array([], dtype='bool')
    maskSPBoundBinary_CEpos = np.array([], dtype='bool')

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
        

    # CE FOCUSED : systems in CE dataframe are based in the masks

        #Helps to select systems which one of the pairs in the binary is a BH before the first CE formation
        maskCEBH1 = ((stellarTypeCE1*maskCE  == 14) ^ (stellarTypeCE2*maskCE  == 14))
        #Helps to select systems which one of the pairs in the binary is a BH after the last CE formation
        maskCEBH2 = ((stellarTypeCE1*mask_unique  == 14) ^ (stellarTypeCE2*mask_unique  == 14)) 
        #Selects systems which both of the pairs in the binary are BH before CE formation
        maskCEBBH1 = ((stellarTypeCE1*maskCE  == 14) & (stellarTypeCE2*maskCE  == 14))
        #Selects systems which both of the pairs in the binary are BH after CE formation
        maskCEBBH2 = ((stellarTypeCE1*mask_unique  == 14) & (stellarTypeCE2*mask_unique  == 14))

    # DC FOCUSED : systems in DC dataframe are based in the masks
        maskDC_CE = np.in1d(seedsDC, seedsCE)
        #DCs in all SPs
        maskDCBH1 = ((stellarTypeDC1  == 14) ^ (stellarTypeDC2  == 14)) 
        maskDCBH1pr = ((stellarTypeDC1  == 14) & (stellarTypeDC2  != 14)) 
        maskDCBH1sc = ((stellarTypeDC1  != 14) & (stellarTypeDC2  == 14)) 

        maskDCBBH1 = ((stellarTypeDC1  == 14) & (stellarTypeDC2  == 14))

        #DCs only CE positive
        maskDCBH1_CEpos = ((stellarTypeDC1*maskDC_CE  == 14) ^ (stellarTypeDC2*maskDC_CE  == 14))
        maskDCBH1pr_CEpos = ((stellarTypeDC1*maskDC_CE  == 14) & (stellarTypeDC2*maskDC_CE  != 14))
        maskDCBH1sc_CEpos= ((stellarTypeDC1*maskDC_CE  != 14) & (stellarTypeDC2*maskDC_CE  == 14))

        maskDCBBH1_CEpos = ((stellarTypeDC1*maskDC_CE == 14) & (stellarTypeDC2*maskDC_CE  == 14))

    # SP FOCUSED : systems in SP dataframe are based in the masks
        mask = np.in1d(seedsSP, seedsCE)

        maskSPBH1 = ((stellarTypeSP1  == 14) ^ (stellarTypeSP2  == 14)) 
        maskSPBH1pr = ((stellarTypeSP1  == 14) & (stellarTypeSP2  != 14)) 
        maskSPBH1sc = ((stellarTypeSP1  != 14) & (stellarTypeSP2  == 14)) 

        maskSPBH2 = ((stellarTypeSP1  == 14) & (stellarTypeSP2  == 14))
        
        #SPs only CE positive
        maskSPBH1_CEpos = ((stellarTypeSP1*mask  == 14) ^ (stellarTypeSP2*mask  == 14)) 
        maskSPBH1pr_CEpos = ((stellarTypeSP1*mask  == 14) & (stellarTypeSP2*mask  != 14)) 
        maskSPBH1sc_CEpos = ((stellarTypeSP1*mask  != 14) & (stellarTypeSP2*mask  == 14)) 

        maskSPBH2_CEpos= ((stellarTypeSP1*mask  == 14) & (stellarTypeSP2*mask  == 14)) 

        maskSPunb = (status == 14) & maskSPBH1
        maskSPdco = (status == 11) & maskSPBH1
        maskSPmrgr = (status == 9)  &  maskSPBH1
        maskSPunb_CEpos = (status*mask == 14) & maskSPBH1_CEpos 
        maskSPdco_CEpos = (status*mask == 11) & maskSPBH1_CEpos 
        maskSPmrgr_CEpos = (status*mask == 9) & maskSPBH1_CEpos 
        #Note: if you were used status*mask then maskSPunb_CEpos would have the same len with CE positive systems, and on the plots you would use mP_ins to mask.
        maskSPBndBnry = (maskSPBH1 == 1) & (maskSPdco == 1) & (maskSPunb == 0)
        maskSPBndBnry_CEpos = (maskSPBH1_CEpos == 1) & (maskSPdco_CEpos == 1) & (maskSPunb_CEpos == 0)
        print("pre CE maskpreCEBH: ", np.sum(maskCEBH1), "post CE maskpstCEBH: ", np.sum(maskCEBH2))


        # maskDCBH1 = ((stellarTypeDC1  == 14) | (stellarTypeDC2  == 14)) 
        # maskDCBH2 = ((stellarTypeDC1  == 14) & (stellarTypeDC2  == 14))
        # maskDCBH1_CEpos = ((stellarTypeDC1*mask  == 14) | (stellarTypeDC2*mask  == 14)) 
        # maskDCBH2_CEpos= ((stellarTypeDC1*mask  == 14) & (stellarTypeDC2*mask  == 14)) 

        # maskDCunb = (status == 14)
        # maskDCdco = (status == 11)
        # maskDCmrgr = (status == 9) 
        # maskDCunb_CEpos = (status*mask == 14) 
        # maskDCdco_CEpos = (status*mask == 11) 
        # maskDCmrgr_CEpos = (status*mask == 9) 


        mP_innm = SP['Mass@ZAMS(1)'][()]
        mS_innm = SP['Mass@ZAMS(2)'][()]
        mP_in = SP['Mass@ZAMS(1)'][()]*mask
        mS_in = SP['Mass@ZAMS(2)'][()]*mask

        mP_preCE = CE["Mass(1)<CE"][()]*maskCE
        mS_preCE = CE["Mass(2)<CE"][()]*maskCE

        mP_pstCE = CE["Mass(1)>CE"][()]*mask_unique
        mS_pstCE = CE["Mass(2)>CE"][()]*mask_unique

        mP_DC = DC['Mass(1)']
        mS_DC = DC['Mass(2)']

        mP_DC_CE = DC['Mass(1)']*maskDC_CE
        mS_DC_CE = DC['Mass(2)']*maskDC_CE


        # mP_in = SP['Mass@ZAMS(1)'][()]
        # mS_in = SP['Mass@ZAMS(2)'][()]

        # mP_preCE = CE["Mass(1)<CE"][()]*maskCE
        # mS_preCE = CE["Mass(2)<CE"][()]*maskCE

        # mP_pstCE = CE["Mass(1)>CE"][()]*mask_unique
        # mS_pstCE = CE["Mass(2)>CE"][()]*mask_unique

        # mtP_in = SP['Mass@ZAMS(1)'][()][maskCEBH1]
        # mtS_in = SP['Mass@ZAMS(2)'][()][maskCEBH1]

        # mtP_preCE = CE["Mass(1)<CE"][()][maskCEBH1]
        # mtS_preCE = CE["Mass(2)<CE"][()][maskCEBH1]

        # mtP_pstCE = CE["Mass(1)>CE"][()][maskBH2]
        # mtS_pstCE = CE["Mass(2)>CE"][()][maskBH2]


        # orbital_period = SP['Orbital_Period'][()]
        semax_innm = SP['SemiMajorAxis@ZAMS'][()]
        semax_in = SP['SemiMajorAxis@ZAMS'][()]*mask
        semax_preCE = CE["SemiMajorAxis<CE"][()]*maskCE
        semax_pstCE = CE["SemiMajorAxis>CE"][()]*mask_unique
        semax_DC = DC['SemiMajorAxis@DCO'][()]
        semax_DC_CE = DC['SemiMajorAxis@DCO'][()]*maskDC_CE

        mP_innms.extend(mP_innm)
        mP_ins.extend(mP_in)
        mP_preCEs.extend(mP_preCE)
        mP_pstCEs.extend(mP_pstCE)
        mP_DCs.extend(mP_DC)
        mP_DCs_CEs.extend(mP_DC_CE)

        mS_innms.extend(mS_innm)
        mS_ins.extend(mS_in)
        mS_preCEs.extend(mS_preCE)
        mS_pstCEs.extend(mS_pstCE)
        mS_DCs.extend(mS_DC)
        mS_DCs_CEs.extend(mS_DC_CE)


        # orbital_periods.extend(orbital_period)
        # orbital_periods_masked.extend(orbital_period*mask)
        semax_innms.extend(semax_innm)
        semax_ins.extend(semax_in)
        semax_DCs.extend(semax_DC)
        semax_DCs_CEs.extend(semax_DC_CE)
        semax_preCEs.extend(semax_preCE)
        semax_pstCEs.extend(semax_pstCE)

        maskSP = np.concatenate([maskSP, mask])
        maskSPCE = np.concatenate([maskSPCE, maskCE])
        final_status = np.concatenate([final_status, status])

        # np.concatenate([maskCEpre, maskCE)
        # np.concatenate([maskCEpst, mask_unique)

        maskpreCEBH = np.concatenate([maskpreCEBH, maskCEBH1])
        maskpstCEBH = np.concatenate([maskpstCEBH, maskCEBH2])
        maskpreCEBBH = np.concatenate([maskpreCEBBH, maskCEBBH1])
        maskpstCEBBH = np.concatenate([maskpstCEBBH, maskCEBBH2])
        maskSPBH = np.concatenate([maskSPBH, maskSPBH1])
        maskSPBHpr = np.concatenate([maskSPBHpr, maskSPBH1pr])
        maskSPBHsc = np.concatenate([maskSPBHsc, maskSPBH1sc])
        maskSPBBH = np.concatenate([maskSPBBH, maskSPBH2])

        maskDCBH = np.concatenate([maskDCBH, maskDCBH1])
        maskDCBHpr = np.concatenate([maskDCBHpr, maskDCBH1pr])
        maskDCBHsc = np.concatenate([maskDCBHsc, maskDCBH1sc])
        maskDCBBH = np.concatenate([maskDCBBH, maskDCBBH1])

        maskSPUnbound = np.concatenate([maskSPUnbound, maskSPunb])
        maskSPDCO = np.concatenate([maskSPDCO, maskSPdco])
        maskSPMerg = np.concatenate([maskSPMerg, maskSPmrgr])

        maskSPBH_CEpos = np.concatenate([maskSPBH_CEpos, maskSPBH1_CEpos])
        maskSPBHpr_CEpos = np.concatenate([maskSPBHpr_CEpos, maskSPBH1pr_CEpos])
        maskSPBHsc_CEpos = np.concatenate([maskSPBHsc_CEpos, maskSPBH1sc_CEpos])
        maskSPBBH_CEpos = np.concatenate([maskSPBBH_CEpos, maskSPBH2_CEpos])

        maskDCBH_CEpos = np.concatenate([maskDCBH_CEpos, maskDCBH1_CEpos])
        maskDCBHpr_CEpos = np.concatenate([maskDCBHpr_CEpos, maskDCBH1pr_CEpos])
        maskDCBHsc_CEpos = np.concatenate([maskDCBHsc_CEpos, maskDCBH1sc_CEpos])
        maskDCBBH_CEpos = np.concatenate([maskDCBBH_CEpos, maskDCBBH1_CEpos])

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

    L = [ "\n Number of systems: ",str(len(maskSP)),  "\n Number of systems undergo CE: ", str(np.sum(maskSP)), "\n Number of systems undergo CE more than once: ",str(np.sum(maskCEMult)),  
        "\n Percentage of systems undergo CE: ",str((np.sum(maskSP)*100)/len(maskSP)),  " \n Percentage of systems undergo CE more than once: ",str((np.sum(maskCEMult)*100/len(maskCEMult))), 
        "\n Number of BH in all Systems: ", str(np.sum(maskSPBH)),  "\n Number of BH in CE positive systems: ", str(np.sum(maskSPBH_CEpos)),  
        "\n Percentage of BH in all systems: ", str(np.sum(maskSPBH)*100/len(maskSPBH)),  "\n Percentage of BH in CE positive systems: ", str(np.sum(maskSPBH_CEpos)*100/len(maskSPBH_CEpos)),
        "\n Number of BH (only primary) in all Systems: ", str(np.sum(maskSPBHpr)),  "\n Number of BH (only primary) in CE positive systems: ", str(np.sum(maskSPBHpr_CEpos)),  
        "\n Percentage of BH (only primary) in all systems: ", str(np.sum(maskSPBHpr)*100/len(maskSPBHpr)),  "\n Percentage of BH (only primary) in CE positive systems: ", str(np.sum(maskSPBHpr_CEpos)*100/len(maskSPBHpr_CEpos)),
        "\n Number of BH (only secondary) in all Systems: ", str(np.sum(maskSPBHsc)),  "\n Number of BH (only secondary) in CE positive systems: ", str(np.sum(maskSPBHsc_CEpos)),  
        "\n Percentage of BH (only secondary) in all systems: ", str(np.sum(maskSPBHsc)*100/len(maskSPBHsc)),  "\n Percentage of BH (only secondary) in CE positive systems: ", str(np.sum(maskSPBHsc_CEpos)*100/len(maskSPBHsc_CEpos)),        
        "\n Number of BBH in all Systems: ", str(np.sum(maskSPBBH)),  "\n Number of BBH CE positive systems: ", str(np.sum(maskSPBBH_CEpos)),  
        "\n Percentage of BBH in all systems: ", str(np.sum(maskSPBBH)*100/len(maskSPBBH)),  "\n Percentage of BBH in CE positive systems: ", str(np.sum(maskSPBBH_CEpos)*100/len(maskSPBBH_CEpos)),  
        "\n Number of unbound binaries in all systems: ", str(np.sum(maskSPUnbound)), "\n Number of unbound binaries in CE positive systems: ", str(np.sum(maskSPUnbound_CEpos)), 
        "\n Percentage of unbound binaries in all systems: ", str(np.sum(maskSPUnbound)*100/len(maskSPUnbound)), "\n Percentage of unbound binaries in CE positive systems: ", str(np.sum(maskSPUnbound_CEpos)*100/len(maskSPUnbound_CEpos)),  
        "\n Number of double compact objects in all systems: ", str(np.sum(maskSPDCO)),  "\n Number of double compact objects in CE positive systems: ", str(np.sum(maskSPDCO_CEpos)), 
        "\n Percentage of double compact objects in all systems: ", str(np.sum(maskSPDCO)*100/len(maskSPDCO)),  "\n Percentage of double compact objects in CE positive systems: ", str(np.sum(maskSPDCO_CEpos)*100/len(maskSPDCO_CEpos)), 
        "\n Number of mergers in all systems: ", str(np.sum(maskSPMerg)),  "\n Number of mergers in CE positive systems: ", str(np.sum(maskSPMerg_CEpos)),  
        "\n Percentage of mergers in all systems: ", str(np.sum(maskSPMerg)*100/len(maskSPMerg)), "\n Percentage of mergers in CE positive systems: ", str(np.sum(maskSPMerg_CEpos)*100/len(maskSPMerg_CEpos)), 
        "\n Number of bound binaries with at least one BH in all systems: ", str(np.sum(maskSPBoundBinary)), "\n Number of bound binaries with at least one BH in CE positive systems: ", str(np.sum(maskSPBoundBinary_CEpos)),
        "\n Percentage of bound binaries with at least one BH in all systems: ", str(np.sum(maskSPBoundBinary)*100/len(maskSPBoundBinary)), "\n Percentage of bound binaries with at least one BH in CE positive systems: ", str(np.sum(maskSPBoundBinary_CEpos*100/len(maskSPBoundBinary_CEpos))),
        
        "\n Number of BH in DCOs: ", str(np.sum(maskDCBH)),  "\n Number of BH in CE positive DCOs: ", str(np.sum(maskDCBH_CEpos)),  
        "\n Percentage of BH in DCOs: ", str(np.sum(maskDCBH)*100/len(maskDCBH)),  "\n Percentage of BH in CE positive DCOs: ", str(np.sum(maskDCBH_CEpos)*100/len(maskDCBH_CEpos)),
        "\n Number of BH (only primary) in all Systems: ", str(np.sum(maskDCBHpr)),  "\n Number of BH (only primary) in CE positive DCOs: ", str(np.sum(maskDCBHpr_CEpos)),  
        "\n Percentage of BH (only primary) in DCOs: ", str(np.sum(maskDCBHpr)*100/len(maskDCBHpr)),  "\n Percentage of BH (only primary) in CE positive DCOs: ", str(np.sum(maskDCBHpr_CEpos)*100/len(maskDCBHpr_CEpos)),
        "\n Number of BH (only secondary) in all Systems: ", str(np.sum(maskDCBHsc)),  "\n Number of BH (only secondary) in CE positive DCOs: ", str(np.sum(maskDCBHsc_CEpos)),  
        "\n Percentage of BH (only secondary) in DCOs: ", str(np.sum(maskDCBHsc)*100/len(maskDCBHsc)),  "\n Percentage of BH (only secondary) in CE positive DCOs: ", str(np.sum(maskDCBHsc_CEpos)*100/len(maskDCBHsc_CEpos)),        
        "\n Number of BBH in DCOs: ", str(np.sum(maskDCBBH)),  "\n Number of BBH CE positive DCOs: ", str(np.sum(maskDCBBH_CEpos)),  
        "\n Percentage of BBH in DCOs: ", str(np.sum(maskDCBBH)*100/len(maskDCBBH)),  "\n Percentage of BBH in CE positive DCOs: ", str(np.sum(maskDCBBH_CEpos)*100/len(maskDCBBH_CEpos)),  
        
        "\n All:"
        "\n Simulation completed: ", str(np.sum(final_status == 1)), 
        "\n Evolution stopped because an error occurred: ",  str(np.sum(final_status == 2)),
        "\n Allowed time exceeded: ",  str(np.sum(final_status == 3)),
        "\n Allowed timesteps exceeded: ",  str(np.sum(final_status == 4)),
        "\n SSE error for one of the constituent stars: ",  str(np.sum(final_status == 5)),
        "\n Error evolving binary: ",  str(np.sum(final_status == 6)),
        "\n Time exceeded DCO merger time: ",  str(np.sum(final_status == 7)),
        "\n Stars touching: ",  str(np.sum(final_status == 8)),
        "\n Stars merged: ",str(np.sum(final_status == 9)),
        "\n Stars merged at birth: ",  str(np.sum(final_status == 10)),
        "\n DCO formed: ",  str(np.sum(final_status == 11)),
        "\n Double White Dwarf formed: ",  str(np.sum(final_status == 12)),
        "\n Massless remnant formed: ",  str(np.sum(final_status == 13)),
        "\n Unbound binary: ", str(np.sum(final_status == 14)),

        "\n Single BH:"
        "\n Simulation completed: ", str(np.sum(final_status*maskSPBH == 1)), 
        "\n Evolution stopped because an error occurred: ",  str(np.sum(final_status*maskSPBH == 2)),
        "\n Allowed time exceeded: ",  str(np.sum(final_status*maskSPBH == 3)),
        "\n Allowed timesteps exceeded: ",  str(np.sum(final_status*maskSPBH == 4)),
        "\n SSE error for one of the constituent stars: ",  str(np.sum(final_status*maskSPBH == 5)),
        "\n Error evolving binary: ",  str(np.sum(final_status*maskSPBH == 6)),
        "\n Time exceeded DCO merger time: ",  str(np.sum(final_status*maskSPBH == 7)),
        "\n Stars touching: ",  str(np.sum(final_status*maskSPBH == 8)),
        "\n Stars merged: ",  str(np.sum(final_status*maskSPBH == 9)),
        "\n Stars merged at birth: ",  str(np.sum(final_status*maskSPBH == 10)),
        "\n DCO formed: ",  str(np.sum(final_status*maskSPBH == 11)),
        "\n Double White Dwarf formed: ",  str(np.sum(final_status*maskSPBH == 12)),
        "\n Massless remnant formed: ",  str(np.sum(final_status*maskSPBH == 13)),
        "\n Unbound binary: ", str(np.sum(final_status*maskSPBH == 14)),

        ]

    f.writelines(L)
    f.close()
    print("maskSPBH: ", maskSPBH, " maskSPBBH: ", maskSPBBH)
    # maskSP = np.invert(maskSP)
    # maskpstCEBH = np.invert(maskpstCEBH)
    # maskpreCEBH = np.invert(maskpreCEBH)
    # maskpreCEBBH = np.invert(maskpreCEBBH)
    # maskpstCEBBH = np.invert(maskpstCEBBH)

    # maskSPBH = np.invert(maskSPBH)
    # maskSPBHpr = np.invert(maskSPBHpr)
    # maskSPBHsc = np.invert(maskSPBHsc)
    # maskSPBBH = np.invert(maskSPBBH)

    # maskDCBH = np.invert(maskDCBH)
    # maskDCBHpr = np.invert(maskDCBHpr)
    # maskDCBHsc = np.invert(maskDCBHsc)
    # maskDCBBH = np.invert(maskDCBBH)

    # maskSPBBH = np.invert(maskSPBBH)
    # maskSPDCO = np.invert(maskSPDCO)
    # maskSPUnbound = np.invert(maskSPUnbound)
    # maskSPMerg = np.invert(maskSPMerg)
    # maskSPBoundBinary = np.invert(maskSPBoundBinary)

    # maskSPBH_CEpos = np.invert(maskSPBH_CEpos)
    # maskSPBHpr_CEpos = np.invert(maskSPBHpr_CEpos)
    # maskSPBHsc_CEpos = np.invert(maskSPBHsc_CEpos)
    # maskSPBBH_CEpos = np.invert(maskSPBBH_CEpos)

    # maskDCBH_CEpos = np.invert(maskDCBH_CEpos)
    # maskDCBHpr_CEpos = np.invert(maskDCBHpr_CEpos)
    # maskDCBHsc_CEpos = np.invert(maskDCBHsc_CEpos)
    # maskDCBBH_CEpos = np.invert(maskDCBBH_CEpos)

    # maskSPDCO_CEpos = np.invert(maskSPDCO_CEpos)
    # maskSPUnbound_CEpos = np.invert(maskSPUnbound_CEpos)
    # maskSPMerg_CEpos = np.invert(maskSPMerg_CEpos)
    # maskSPBoundBinary_CEpos = np.invert(maskSPBoundBinary_CEpos)

    # print("mP_ins BH: ", len(ma.masked_array(mP_ins, mask=maskSPBH_CEpos)), " mP_ins BBH : ", len(ma.masked_array(mP_ins, mask=maskSPBBH_CEpos)))

##UNCOMMENT
    # arrays_to_save = [mP_innms, mP_ins, mP_preCEs, mP_pstCEs, mP_DCs, mP_DCs_CEs, mS_innms, mS_ins, mS_preCEs, mS_pstCEs, mS_DCs, mS_DCs_CEs, semax_innms, semax_ins, semax_preCEs, semax_pstCEs, semax_DCs, semax_DCs_CEs,
    #                   maskpreCEBH, maskpstCEBH, maskpreCEBBH, maskpstCEBBH, maskSPBH, maskSPBHpr, maskSPBHsc, maskSPBBH, maskDCBH, maskDCBHpr, maskDCBHsc, maskDCBBH, 
    #                   maskSPUnbound, maskSPDCO, maskSPMerg, 
    #                   maskSPBH_CEpos, maskSPBHpr_CEpos, maskSPBHsc_CEpos, maskSPBBH_CEpos, maskDCBH_CEpos, maskDCBHpr_CEpos, maskDCBHsc_CEpos, maskDCBBH_CEpos, maskSPBBH_CEpos, maskSPUnbound_CEpos, maskSPDCO_CEpos, maskSPMerg_CEpos, maskCEMult, maskSPBoundBinary, maskSPBoundBinary_CEpos, final_status]

    # filenames =  ['mP_innms', 'mP_ins', 'mP_preCEs', 'mP_pstCEs', 'mP_DCs', 'mP_DCs_CEs','mS_innms', 'mS_ins', 'mS_preCEs', 'mS_pstCEs', 'mS_DCs', 'mS_DCs_CEs', 'semax_innms', 'semax_ins', 'semax_preCEs', 'semax_pstCEs', 'semax_DCs','semax_DCs_CEs', 
    #               'maskpreCEBH', 'maskpstCEBH', 'maskpreCEBBH', 'maskpstCEBBH', 'maskSPBH', 'maskSPBHpr', 'maskSPBHsc', 'maskSPBBH', 'maskDCBH', 'maskDCBHpr', 'maskDCBHsc', 'maskDCBBH',
    #               'maskSPUnbound', 'maskSPDCO', 'maskSPMerg', 
    #               'maskSPBH_CEpos', 'maskSPBHpr_CEpos', 'maskSPBHsc_CEpos', 'maskSPBBH_CEpos', 'maskDCBH_CEpos', 'maskDCBHpr_CEpos', 'maskDCBHsc_CEpos', 'maskDCBBH_CEpos', 'maskSPBBH_CEpos', 'maskSPUnbound_CEpos', 'maskSPDCO_CEpos', 'maskSPMerg_CEpos', 'maskCEMult', 'maskSPBoundBinary', 'maskSPBoundBinary_CEpos', 'finalStatus']

    # i = 0
    # for f in arrays_to_save:
    #     name = filenames[i]
    #     np.savetxt("/data/a.saricaoglu/Files/" + model + "/" +  str(c.strftime("%m.%d")) + "/" + mode + "_" + name + "_" + str(c.strftime("%m.%d")) + ".txt", f)
    #     i+=1
###UNCOMMENT
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
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model + "_hist_" + current_time + ".png")
    # plt.close()

    # x_scale = np.logspace(np.log10(1.0), np.log10(max(mP_ins)), 50)
    # y_scale = np.logspace(np.log10(1.0), np.log10(max(mS_ins)), 50)

    # plt.subplot(3,3,1)
    # plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + "\n Evolution Status vs Masses in CE Positive Systems")
    # plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPBH_CEpos)-ma.masked_array(mP_ins, mask= maskSPBBH_CEpos)), 
    #            abs(ma.masked_array(mS_ins, mask=maskSPBH_CEpos)-ma.masked_array(mS_ins, mask=maskSPBBH_CEpos)),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial  [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "(BH-BBH) Difference")

    # plt.subplot(3,3,2)
    # plt.hist2d(abs(ma.masked_array(mP_preCEs,mask=maskpreCEBH)-ma.masked_array(mP_preCEs, mask=maskpreCEBBH)), 
    #            abs(ma.masked_array(mS_preCEs, mask=maskpreCEBH)-ma.masked_array(mS_preCEs, mask=maskpreCEBBH)),bins=(x_scale, y_scale),)
    # plt.xlabel("$M_p$ before CE [$M_{sun}$]")
    # plt.ylabel("$M_s$ before CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "(BH-BBH) Difference")

    # plt.subplot(3,3,3)
    # plt.hist2d(abs(ma.masked_array(mP_pstCEs, mask=maskpstCEBH)-ma.masked_array(mP_pstCEs, mask=maskpstCEBBH)),
    #             (ma.masked_array(mS_pstCEs, mask=maskpstCEBH)-ma.masked_array(mS_pstCEs, mask=maskpstCEBBH)),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ after CE [$M_{sun}$]")
    # plt.ylabel("$M_s$ after CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "(BH-BBH) Difference")

    # plt.subplot(3,3,4)
    # plt.hist2d(ma.masked_array(mP_ins, mask= maskSPBH_CEpos), ma.masked_array(mS_ins, mask=maskSPBH_CEpos),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "Black Hole occurance")

    # plt.subplot(3,3,5)
    # plt.hist2d(ma.masked_array(mP_ins, mask= maskSPBBH_CEpos), ma.masked_array(mS_ins, mask=maskSPBBH_CEpos),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial CE [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "Binary Black Hole occurance")

    # plt.subplot(3,3,6)
    # plt.hist2d(ma.masked_array(mP_ins, mask=maskSPUnbound_CEpos), ma.masked_array(mS_ins, mask=maskSPUnbound_CEpos),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "Unbound Binary (1 BH) occurance")

    # plt.subplot(3,3,7)
    # plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPDCO_CEpos)-ma.masked_array(mP_ins, mask= maskSPBBH_CEpos)),
    #             abs(ma.masked_array(mS_ins, mask=maskSPDCO_CEpos)-ma.masked_array(mS_ins, mask= maskSPBBH_CEpos)),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "(DCO-BBH) Difference")

    # plt.subplot(3,3,8)
    # plt.hist2d(abs(ma.masked_array(mP_ins, mask=maskSPUnbound_CEpos)-ma.masked_array(mP_ins, mask=maskSPMerg_CEpos)), 
    #            abs(ma.masked_array(mS_ins, mask=maskSPUnbound_CEpos)- ma.masked_array(mS_ins, mask=maskSPMerg_CEpos)),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial[$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "(UBB-Merg) Difference")

    # plt.subplot(3,3,9)
    # plt.hist2d(ma.masked_array(mP_ins, mask=maskSPDCO_CEpos), ma.masked_array(mS_ins, mask=maskSPDCO_CEpos),bins=(x_scale, y_scale))
    # plt.xlabel("$M_p$ initial [$M_{sun}$]")
    # plt.ylabel("$M_s$ initial [$M_{sun}$]")
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.colorbar(label = "DCO (except BBH) occurance")


    # plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.45)
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model + "_2dhist_" + current_time + ".png")
    # plt.close()

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
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model + "_hist_" + current_time + ".png")
    # plt.close()

    # y_scale2 = np.logspace(np.log10(1.0), np.log10(max(semax_innms)), 50)



    # plt.subplot(3,2,1)
    # plt.suptitle(model+ " IMF, Population size: " + str(len((mP_innms))) + "\n Evolution Status vs Semimajor Axis in CE Positive Systems")
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPBH_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" Single BH Binaries ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplot(3,2,2)
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPBBH_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" BH Binaries ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplot(3,2,3)
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPBoundBinary_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" Bound Single BH Binaries ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplot(3,2,4)
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPDCO_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" DCOs w/ Single BH ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplot(3,2,5)
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPUnbound_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" Unbound Sys. w/ Single BH ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplot(3,2,6)
    # plt.hist([ma.masked_array(semax_ins, mask=maskSPMerg_CEpos)],bins=50)
    # plt.xlabel("$a$ initial [AU]")
    # plt.ylabel(" Mergers w/ Single BH ")
    # plt.yscale("log")
    # #plt.xscale("log")

    # plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.32, wspace=0.45)
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model + "_histsemax_" + current_time + ".png")
    # plt.close()
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
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model +"_scatter_" + current_time + ".png")
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
    # plt.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/Single_BH/" + model + "_boundbh_" + current_time + ".png")
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

