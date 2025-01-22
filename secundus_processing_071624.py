#### PARENT SCRIPT: kroupa_general_multplmask_final_SPsDCOs.py
#----> a refined version for Github push. (16.07.24)
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
print("start time :", current_time)

models = ["Kroupa"]

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

    # Keeps corresponding numerical values 

    SPs = []
    DCs = []

    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []       
    STELLARTYPEDC1   =  []
    STELLARTYPEDC2   =  [] 

    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    MASSDC1 = []
    MASSDC2 = []

    SEMIMAJORAXISSP = []
    SEMIMAJORAXISDC = []

    # Boolen arrays for masking.
    MASKSPBH1 = np.array([], dtype='bool')
    MASKSPBH2 = np.array([], dtype='bool')
    MASKDCBH1 = np.array([], dtype='bool')
    MASKDCBH2 = np.array([], dtype='bool')

    MASKSPunb = np.array([], dtype='bool')
    MASKSPdco = np.array([], dtype='bool')
    MASKSPmrgr = np.array([], dtype='bool')

    MASKDCinSP = np.array([], dtype='bool')
    MASKDCBHNS1 = np.array([], dtype='bool')
    MASKDCBHNS2 = np.array([], dtype='bool')

    MASKDCnonBH =  np.array([], dtype='bool')
    MASKDCSEMAJ = np.array([], dtype='bool')


    f = open("/data/a.saricaoglu/Files/Kroupa/" + model + "_" + str(c.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    f.writelines(["\n","\n Run :", current_time])

    for Data in data_outputs:


        SP = Data['BSE_System_Parameters']
        DC = Data['BSE_Double_Compact_Objects']

        seedsSP = SP['SEED'][()]
        seedsDC = DC['SEED'][()]
        statusSP = SP['Evolution_Status'][()]

        stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
        stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]        
        stellarTypeDC1   =  DC['Stellar_Type(1)'][()]
        stellarTypeDC2   =  DC['Stellar_Type(2)'][()]  

        massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
        massZamsSP2 = SP['Mass@ZAMS(2)'][()] 
        massDCO1 = DC['Mass(1)'][()] 
        massDCO2 = DC['Mass(2)'][()] 

        semimajorAxisSP = SP['SemiMajorAxis@ZAMS'][()] 
        semimajorAxisDC = DC['SemiMajorAxis@DCO'][()] 
## may add eccentricity later on

        maskSPBH1 = ((stellarTypeSP1  == 14) & (stellarTypeSP2  != 14)) 
        maskSPBH2 = ((stellarTypeSP1  != 14) & (stellarTypeSP2  == 14))
        maskDCBH1 = ((stellarTypeDC1  == 14) & (stellarTypeDC2  != 14)) 
        maskDCBH2 = ((stellarTypeDC1  != 14) & (stellarTypeDC2  == 14))       

        maskSPunb = (statusSP == 17)
        maskSPdco = (statusSP == 14)
        maskSPmrgr = (statusSP == 12) 

        maskDCinSP =  np.in1d(seedsSP, seedsDC)

        maskDCnonBH = ((stellarTypeDC1 != 14) & (stellarTypeDC2 != 14)) #none is bh
        maskDCBHNS1 = ((stellarTypeDC1 == 14) & (stellarTypeDC2 == 13)) #primary is bh secondary is ns
        maskDCBHNS2 = ((stellarTypeDC1 == 13) & (stellarTypeDC2 == 14)) #primary is ns secondary is bh

        maskDCsemaj = (semimajorAxisDC <= 0.5)

        SPs.extend(seedsSP)
        DCs.extend(seedsDC)

        STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
        STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
        STELLARTYPESP1.extend(stellarTypeSP1)
        STELLARTYPESP2.extend(stellarTypeSP2)       
        STELLARTYPEDC1.extend(stellarTypeDC1)
        STELLARTYPEDC2.extend(stellarTypeDC2) 

        MASSZAMSSP1.extend(massZamsSP1)
        MASSZAMSSP2.extend(massZamsSP2)
        MASSDC1.extend(massDCO1)
        MASSDC2.extend(massDCO2)

        SEMIMAJORAXISSP.extend(semimajorAxisSP)
        SEMIMAJORAXISDC.extend(semimajorAxisDC)

        MASKSPBH1 = np.concatenate([MASKSPBH1, maskSPBH1])
        MASKSPBH2 = np.concatenate([MASKSPBH2, maskSPBH2])
        MASKDCBH1 =np.concatenate([MASKDCBH1, maskDCBH1])
        MASKDCBH2 = np.concatenate([MASKDCBH2, maskDCBH2])

        MASKSPunb = np.concatenate([MASKSPunb, maskSPunb])
        MASKSPdco = np.concatenate([MASKSPdco, maskSPdco])
        MASKSPmrgr = np.concatenate([MASKSPmrgr, maskSPmrgr])


        MASKDCinSP = np.concatenate([MASKDCinSP, maskDCinSP])
        MASKDCBHNS1 = np.concatenate([MASKDCBHNS1, maskDCBHNS1])
        MASKDCBHNS2 = np.concatenate([MASKDCBHNS2, maskDCBHNS2])

        MASKDCnonBH = np.concatenate([MASKDCnonBH, maskDCnonBH])
        MASKDCSEMAJ =  np.concatenate([MASKDCSEMAJ, maskDCsemaj])
      

        Data.close()

    L = ["\n Number of systems: ",str(len(MASSZAMSSP1)),  "\n Number of systems forming DCO (via SP status): ", str(np.sum(MASKSPdco)), "\n Number of systems forming DCO (via DCO mass): ", str(len(MASSDC1)), 
        "\n Number of primary mass BHs (SP): ", str(np.sum(MASKSPBH1)), "\n Number of primary mass BHs (DC): ", str(np.sum(MASKDCBH1)),  
        "\n Number of secondary mass BHs (SP): ", str(np.sum(MASKSPBH2)), "\n Number of secondary mass BHs (DC): ", str(np.sum(MASKDCBH2)),     
        "\n Number of unbound binaries (SP): ", str(np.sum(MASKSPunb)), 
        "\n Number of double compact objects (SP): ", str(np.sum(MASKSPdco)),  
        "\n Number of mergers (SP): ", str(np.sum(MASKSPmrgr)),
        "\n Number of BH-NS binaries (DC): ", str(np.sum(MASKDCBHNS1)), "\n Number of NS-BH binaries (DC): ", str(np.sum(MASKDCBHNS2)),
        "\n Number of systems without any BHs (DC): ", str(np.sum(MASKDCnonBH)),
        "\n Number of systems with semimajor axis < 0.5 AU (DC): ", str(np.sum(MASKDCSEMAJ))
        ]

    f.writelines(L)
    f.close()

    arrays_to_save = [SPs,DCs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEDC1,STELLARTYPEDC2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSDC1,MASSDC2,SEMIMAJORAXISSP,SEMIMAJORAXISDC,MASKSPBH1,MASKSPBH2,MASKDCBH1,MASKDCBH2,
                    MASKSPunb,MASKSPdco,MASKSPmrgr,MASKDCinSP,MASKDCBHNS1,MASKDCBHNS2, MASKDCnonBH, MASKDCSEMAJ]
    filenames = ["SPs","DCs","STELLARTYPEZAMSSP1","STELLARTYPEZAMSSP2","STELLARTYPESP1","STELLARTYPESP2","STELLARTYPEDC1","STELLARTYPEDC2",
                "MASSZAMSSP1","MASSZAMSSP2","MASSDC1","MASSDC2","SEMIMAJORAXISSP","SEMIMAJORAXISDC","MASKSPBH1","MASKSPBH2","MASKDCBH1","MASKDCBH2",
                "MASKSPunb","MASKSPdco","MASKSPmrgr","MASKDCinSP","MASKDCBHNS1","MASKDCBHNS2", "MASKDCnonBH", "MASKDCSEMAJ"]
    i = 0
    if not os.path.exists("/data/a.saricaoglu/Files/" + model + "/" +  str(c.strftime("%m.%d"))): 
        os.makedirs("/data/a.saricaoglu/Files/" + model + "/" +  str(c.strftime("%m.%d"))) 
    for f in arrays_to_save:
        name = filenames[i]
        np.savetxt("/data/a.saricaoglu/Files/" + model + "/" +  str(c.strftime("%m.%d")) + "/" + name + "_" + str(c.strftime("%m.%d")) + ".txt", f)
        i+=1

c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("finish time :", current_time)