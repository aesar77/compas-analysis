#### PARENT SCRIPT: secundus_processng_081424.py
#----> this version captures the system just before and after common envelope (CE) event.
import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time as t                     # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  # for plotting
import matplotlib
import calculations as calc      # functions from calculations.py

matplotlib.use("Agg")

# Displays Time
s = datetime.now()
starttime = t.ctime(t.time())
start = t.process_time()
start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
print("Start time :", start_time)

# Choose the IMF to process
IMFs = ["Kroupa"]

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with
for IMF in IMFs:
    matplotlib.rcParams['figure.figsize'] = (15,10)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 14
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = '/data/a.saricaoglu/Runs/' + IMF + '/'
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

    # Keeps corresponding numerical values. 

    SPs = []
    CEs = []

    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []       
    STELLARTYPEPRECE1   =  []
    STELLARTYPEPRECE2   =  [] 
    STELLARTYPEPSTCE1   =  []
    STELLARTYPEPSTCE2   =  [] 

    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    MASSPRECE1 = []
    MASSPRECE2 = []
    MASSPSTCE1 = []
    MASSPSTCE2 = []

    SEMIMAJORAXISSP = []
    SEMIMAJORAXISPRECE = []
    SEMIMAJORAXISPSTCE = []

    ORBITALPERIODSP = []
    ORBITALPERIODPRECE = []
    ORBITALPERIODPSTCE= []

    RADIUSPRECE1 = []
    RADIUSPRECE2 = []
    RADIUSPSTCE1 = []
    RADIUSPSTCE2 = []

    COSIPRECE = []
    COSIPSTCE = []

    TIMECE = []

    

    # Boolen arrays for masking.
    MASKSPBH1 = np.array([], dtype='bool')
    MASKSPBH2 = np.array([], dtype='bool')

    MASKPRECE = np.array([], dtype='bool')
    MASKPSTCE = np.array([], dtype='bool')

    MASKSPunb = np.array([], dtype='bool')
    MASKSPdco = np.array([], dtype='bool')
    MASKSPmrgr = np.array([], dtype='bool')

    MASKPRECEBH1 = np.array([], dtype='bool')
    MASKPSTCEBH1 = np.array([], dtype='bool')

    MASKPRECEBH2 = np.array([], dtype='bool')
    MASKPSTCEBH2 = np.array([], dtype='bool')

    MASKPRECEinSP = np.array([], dtype='bool')
    MASKPSTCEinSP = np.array([], dtype='bool')

    MASKPRECEBHNS1 = np.array([], dtype='bool')
    MASKPSTCEBHNS1 = np.array([], dtype='bool')

    MASKPRECEBHNS2 = np.array([], dtype='bool')
    MASKPSTCEBHNS2 = np.array([], dtype='bool')

    MASKPRECEnonBH = np.array([], dtype='bool')
    MASKPSTCEnonBH = np.array([], dtype='bool')

    MASKPRECESEMAJ = np.array([], dtype='bool')
    MASKPSTCESEMAJ = np.array([], dtype='bool')

    MASKPRECERBPER = np.array([], dtype='bool')
    MASKPSTCERBPER = np.array([], dtype='bool')

    MASKPRECEnegCE = np.array([], dtype='bool')
    MASKPSTCEnegCE = np.array([], dtype='bool')

    MASKPRECEposCE = np.array([], dtype='bool')
    MASKPSTCEposCE = np.array([], dtype='bool')

    COSIPRECE = np.array([], dtype='bool')
    COSIPSTCE = np.array([], dtype='bool')

    MASKFIRSTCE = np.array([], dtype='bool')
    MASKLASTCE = np.array([], dtype='bool')

    MASKPRECESEARCHABILITY = np.array([], dtype='bool')
    MASKPSTCESEARCHABILITY = np.array([], dtype='bool')

    MASKENVELOPEJECTION = np.array([], dtype='bool')

    f = open("/data/a.saricaoglu/Files/Kroupa/" + IMF + "_CE_" + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    f.writelines(["\n","\n Run :", start_time])

    i = 1

    for Data in data_outputs:
        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " start time :", current_time)


        SP = Data['BSE_System_Parameters']
        DC = Data['BSE_Double_Compact_Objects']
        CE = Data['BSE_Common_Envelopes']

        seedsSP = SP['SEED'][()]
        seedsDC = DC['SEED'][()]
        seedsCE = CE['SEED'][()]

        eventsCE = CE['CE_Event_Counter'][()]

        statusSP = SP['Evolution_Status'][()]

        stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
        stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]        
        stellarTypepreCE1   =  CE['Stellar_Type(1)<CE'][()]
        stellarTypepreCE2   =  CE['Stellar_Type(2)<CE'][()]  
        stellarTypepstCE1   =  CE['Stellar_Type(1)'][()]
        stellarTypepstCE2   =  CE['Stellar_Type(2)'][()] 
        

        massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
        massZamsSP2 = SP['Mass@ZAMS(2)'][()] 

        # PreCE refers to the mass just before the first CE event, while pstCE refers to the mass just after the last CE event.
        masspreCE1 = CE['Mass(1)<CE'][()] 
        masspreCE2 = CE['Mass(2)<CE'][()] 
        masspstCE1 = CE['Mass(1)>CE'][()] 
        masspstCE2 = CE['Mass(2)>CE'][()]

        semimajorAxisSP = SP['SemiMajorAxis@ZAMS'][()] 
        semimajorAxispreCE = CE['SemiMajorAxis<CE'][()] 
        semimajorAxispstCE = CE['SemiMajorAxis>CE'][()]

        timeCE = CE['Time'][()]
        
        radiuspreCE1 = CE['Radius(1)<CE'][()]
        radiuspstCE1 = CE['Radius(1)>CE'][()] 
        radiuspreCE2 = CE['Radius(2)<CE'][()]
        radiuspstCE2 = CE['Radius(2)>CE'][()]   

        # Chooses first CE events.
        maskfirstCE = (eventsCE == 1)
        # Chooses last CE events
        masklastCE = np.zeros(len(CE['SEED'][()]), dtype='bool')
        last = 0
        for k in range(0, len(CE['SEED'][()])):
            if (CE['SEED'][()][k] == last):
                masklastCE[k-1] = False
                masklastCE[k] = True
            else:
                masklastCE[k] = True
                last = CE['SEED'][()][k]


        orbitalPeriodSP = calc.orbital_period(massZamsSP1, massZamsSP2, semimajorAxisSP)
        orbitalPeriodpreCE = calc.orbital_period(masspreCE1, masspreCE2, semimajorAxispreCE)
        orbitalPeriodpstCE = calc.orbital_period(masspstCE1, masspstCE2, semimajorAxispstCE)
        print(np.shape(orbitalPeriodpreCE))
        cosipreCE = calc.orbital_inclination(radiuspreCE1, radiuspreCE2, semimajorAxispreCE)
        cosipstCE = calc.orbital_inclination(radiuspstCE1, radiuspstCE2, semimajorAxispstCE)
        print(np.shape(cosipreCE))
        maskpreCEsearchability = calc.searchability(cosipreCE)
        maskpstCEsearchability = calc.searchability(cosipstCE)
        print(np.shape(maskpreCEsearchability))
        maskEnvelopejection = calc.envelope_ejection(radiuspstCE1, radiuspstCE2, semimajorAxispstCE)
        print(np.shape(maskEnvelopejection))
        # i = np.random.uniform(0,180,len(seedsCE))



## may add eccentricity later on

        maskSPBH1 = ((stellarTypeSP1  == 14) & (stellarTypeSP2  != 14)) 
        maskSPBH2 = ((stellarTypeSP1  != 14) & (stellarTypeSP2  == 14))
        maskpreCEBH1 = ((stellarTypepreCE1  == 14) & (stellarTypepreCE2  != 14)) 
        maskpreCEBH2 = ((stellarTypepreCE1  != 14) & (stellarTypepreCE2  == 14))       
        maskpstCEBH1 = ((stellarTypepstCE1  == 14) & (stellarTypepstCE2  != 14)) 
        maskpstCEBH2 = ((stellarTypepstCE1  != 14) & (stellarTypepstCE2  == 14))

        maskSPunb = (statusSP == 14)
        maskSPdco = (statusSP == 11)
        maskSPmrgr = (statusSP == 9) 

        maskCEinSP =  np.in1d(seedsSP, seedsCE)

        maskpreCEnonBH = ((stellarTypepreCE1 != 14) & (stellarTypepreCE2 != 14)) #none is bh
        maskpreCEBHNS1 = ((stellarTypepreCE1 == 14) & (stellarTypepreCE2 == 13)) #primary is bh secondary is ns
        maskpreCEBHNS2 = ((stellarTypepreCE1 == 13) & (stellarTypepreCE2 == 14)) #primary is ns secondary is bh

        maskpstCEnonBH = ((stellarTypepstCE1 != 14) & (stellarTypepstCE2 != 14)) #none is bh
        maskpstCEBHNS1 = ((stellarTypepstCE1 == 14) & (stellarTypepstCE2 == 13)) #primary is bh secondary is ns
        maskpstCEBHNS2 = ((stellarTypepstCE1 == 13) & (stellarTypepstCE2 == 14)) #primary is ns secondary is bh

        maskpreCEsemaj = (semimajorAxispreCE <= 3.0)
        maskpstCEsemaj = (semimajorAxispstCE <= 3.0)

        maskpreCEorbper = (orbitalPeriodpreCE <= 30)
        maskpstCEorbper = (orbitalPeriodpstCE <= 30)
        
        maskpreCEnegCE =  np.in1d(seedsCE*maskfirstCE, seedsDC, invert=True)
        maskpreCEposCE =  np.in1d(seedsCE*maskfirstCE, seedsDC)
        maskpstCEnegCE =  np.in1d(seedsCE*masklastCE, seedsDC, invert=True)
        maskpstCEposCE =  np.in1d(seedsCE*masklastCE, seedsDC)

        SPs.extend(seedsSP)
        CEs.extend(seedsCE)

        STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
        STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
        STELLARTYPESP1.extend(stellarTypeSP1)
        STELLARTYPESP2.extend(stellarTypeSP2)       
        STELLARTYPEPRECE1.extend(stellarTypepreCE1)
        STELLARTYPEPRECE2.extend(stellarTypepreCE2) 
        STELLARTYPEPSTCE1.extend(stellarTypepstCE1)
        STELLARTYPEPSTCE2.extend(stellarTypepstCE2) 

        MASSZAMSSP1.extend(massZamsSP1)
        MASSZAMSSP2.extend(massZamsSP2)
        MASSPRECE1.extend(masspreCE1)
        MASSPRECE2.extend(masspreCE2)
        MASSPSTCE1.extend(masspstCE1)
        MASSPSTCE2.extend(masspstCE2)

        SEMIMAJORAXISSP.extend(semimajorAxisSP)
        SEMIMAJORAXISPRECE.extend(semimajorAxispreCE)
        SEMIMAJORAXISPSTCE.extend(semimajorAxispstCE)

        ORBITALPERIODSP.extend(orbitalPeriodSP)
        ORBITALPERIODPRECE.extend(orbitalPeriodpreCE)
        ORBITALPERIODPSTCE.extend(orbitalPeriodpstCE)

        MASKSPunb = np.concatenate([MASKSPunb, maskSPunb])
        MASKSPdco = np.concatenate([MASKSPdco, maskSPdco])
        MASKSPmrgr = np.concatenate([MASKSPmrgr, maskSPmrgr])
        MASKCEinSP = np.concatenate([MASKPRECEinSP, maskCEinSP])


        MASKSPBH1 = np.concatenate([MASKSPBH1, maskSPBH1])
        MASKSPBH2 = np.concatenate([MASKSPBH2, maskSPBH2])
        MASKPRECEBH1 =np.concatenate([MASKPRECEBH1, maskpreCEBH1])
        MASKPRECEBH2 = np.concatenate([MASKPRECEBH2, maskpreCEBH2])

        MASKPRECEBHNS1 = np.concatenate([MASKPRECEBHNS1, maskpreCEBHNS1])
        MASKPRECEBHNS2 = np.concatenate([MASKPRECEBHNS2, maskpreCEBHNS2])

        MASKPRECEnonBH = np.concatenate([MASKPRECEnonBH, maskpreCEnonBH])
        MASKPRECESEMAJ =  np.concatenate([MASKPRECESEMAJ, maskpreCEsemaj])
        MASKPRECERBPER = np.concatenate([MASKPRECERBPER, maskpreCEorbper])
        COSIPRECE = np.concatenate([COSIPRECE, cosipreCE])
        MASKPRECEnegCE =  np.concatenate([MASKPRECEnegCE, maskpreCEnegCE])
        MASKPRECEposCE = np.concatenate([MASKPRECEposCE, maskpreCEposCE])

        MASKPSTCEBH1 =np.concatenate([MASKPSTCEBH1, maskpstCEBH1])
        MASKPSTCEBH2 = np.concatenate([MASKPSTCEBH2, maskpstCEBH2])

        MASKPSTCEBHNS1 = np.concatenate([MASKPSTCEBHNS1, maskpstCEBHNS1])
        MASKPSTCEBHNS2 = np.concatenate([MASKPSTCEBHNS2, maskpstCEBHNS2])

        MASKFIRSTCE = np.concatenate([MASKFIRSTCE, maskfirstCE])
        MASKLASTCE = np.concatenate([MASKLASTCE, masklastCE])

        MASKPSTCEnonBH = np.concatenate([MASKPSTCEnonBH, maskpstCEnonBH])
        MASKPSTCESEMAJ =  np.concatenate([MASKPSTCESEMAJ, maskpstCEsemaj])
        MASKPSTCERBPER = np.concatenate([MASKPSTCERBPER, maskpstCEorbper])
        COSIPSTCE = np.concatenate([COSIPSTCE, cosipstCE])
        MASKPSTCEnegCE =  np.concatenate([MASKPSTCEnegCE, maskpstCEnegCE])
        MASKPSTCEposCE = np.concatenate([MASKPSTCEposCE, maskpstCEposCE])

        MASKPRECESEARCHABILITY = np.concatenate([MASKPRECESEARCHABILITY, maskpreCEsearchability])
        MASKPSTCESEARCHABILITY = np.concatenate([MASKPSTCESEARCHABILITY, maskpstCEsearchability])

        MASKENVELOPEJECTION = np.concatenate([MASKENVELOPEJECTION, maskEnvelopejection])

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " end time :", current_time)
        i += 1

    Data.close()

    L = ["\n Number of systems: ",str(len(MASSZAMSSP1)),  "\n Number of systems forming CE (via SP status): ", str(np.sum(MASKFIRSTCE)), "\n Number of systems forming CE (via CE mass): ", str(len(MASSPRECE1)), 
        "\n Number of primary mass BHs (SP): ", str(np.sum(MASKSPBH1)),   
        "\n Number of secondary mass BHs (SP): ", str(np.sum(MASKSPBH2)),    
        "\n Number of unbound binaries (SP): ", str(np.sum(MASKSPunb)), 
        "\n Number of double compact objects (SP): ", str(np.sum(MASKSPdco)),  
        "\n Number of mergers (SP): ", str(np.sum(MASKSPmrgr)),
        "\n Number of primary mass BHs (pre CE): ", str(np.sum(MASKPRECEBH1)),
         "\n Number of secondary mass BHs (pre CE): ", str(np.sum(MASKPRECEBH2)), 
        "\n Number of BH-NS binaries (pre CE): ", str(np.sum(MASKPRECEBHNS1)), "\n Number of NS-BH binaries (pre CE): ", str(np.sum(MASKPRECEBHNS2)),
        "\n Number of systems without any BHs (pre CE): ", str(np.sum(MASKPRECEnonBH)),
        "\n Number of systems with semimajor axis < 3.0 AU (pre CE): ", str(np.sum(MASKPRECESEMAJ)),
        "\n Number of systems forming DCs (pre CE): ", str(np.sum(MASKPRECEposCE)),
        "\n Number of systems not forming DCs (pre CE): ", str(np.sum(MASKPRECEnegCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND forming DCs (pre CE): ", str(np.sum(MASKPRECESEMAJ*MASKPRECEposCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND NOT forming DCs (pre CE): ", str(np.sum(MASKPRECESEMAJ*MASKPRECEnegCE)),
        "\n Number of systems with orbital period < 30 days AND forming DCs (pre CE): ", str(np.sum(MASKPRECERBPER*MASKPRECEposCE)),
        "\n Number of systems with orbital period < 30 days AND NOT forming DCs (pre CE): ", str(np.sum(MASKPRECERBPER*MASKPRECEnegCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND forming DCs (pre CE): ", str(np.sum(MASKPRECESEARCHABILITY*MASKPRECEposCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND NOT forming DCs (pre CE): ", str(np.sum(MASKPRECESEARCHABILITY*MASKPRECEnegCE)),
         "\n Number of primary mass BHs (pst CE): ", str(np.sum(MASKPSTCEBH1)),
         "\n Number of secondary mass BHs (pst CE): ", str(np.sum(MASKPSTCEBH2)), 
        "\n Number of BH-NS binaries (pst CE): ", str(np.sum(MASKPSTCEBHNS1)), "\n Number of NS-BH binaries (pst CE): ", str(np.sum(MASKPSTCEBHNS2)),
        "\n Number of systems without any BHs (pst CE): ", str(np.sum(MASKPSTCEnonBH)),
        "\n Number of systems with semimajor axis < 3.0 AU (pst CE): ", str(np.sum(MASKPSTCESEMAJ)),
        "\n Number of systems forming DCs (pst CE): ", str(np.sum(MASKPSTCEposCE)),
        "\n Number of systems not forming DCs (pst CE): ", str(np.sum(MASKPSTCEnegCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND forming DCs (pst CE): ", str(np.sum(MASKPSTCESEMAJ*MASKPSTCEposCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND NOT forming DCs (pst CE): ", str(np.sum(MASKPSTCESEMAJ*MASKPSTCEnegCE)),
        "\n Number of systems with orbital period < 30 days AND forming DCs (pst CE): ", str(np.sum(MASKPSTCERBPER*MASKPSTCEposCE)),
        "\n Number of systems with orbital period < 30 days AND NOT forming DCs (pst CE): ", str(np.sum(MASKPSTCERBPER*MASKPSTCEnegCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND forming DCs (pst CE): ", str(np.sum(MASKPSTCESEARCHABILITY*MASKPSTCEposCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND NOT forming DCs (pst CE): ", str(np.sum(MASKPSTCESEARCHABILITY*MASKPSTCEnegCE)),
        "\n Total CE events: ", str(np.sum(eventsCE)),
        "\n Number of systems satisfying envelope ejection condition after the last CE : ", str(np.sum(MASKENVELOPEJECTION)),
        "\n Number of primary HeWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE1 == 10 ])),
        "\n Number of primary COWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE1 == 11 ])),
        "\n Number of primary ONeWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE1 == 12 ])),
        "\n Number of primary HeWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE1 == 10 ])),
        "\n Number of primary COWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE1 == 11 ])),
        "\n Number of primary ONeWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE1 == 12 ])),
        "\n Number of secondary HeWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE2 == 10 ])),
        "\n Number of secondary COWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE2 == 11 ])),
        "\n Number of secondary ONeWD (pre CE) : ", str(np.sum([STELLARTYPEPRECE2 == 12 ])),
        "\n Number of secondary HeWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE2 == 10 ])),
        "\n Number of secondary COWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE2 == 11 ])),
        "\n Number of secondary ONeWD (pst CE) : ", str(np.sum([STELLARTYPEPSTCE2 == 12 ])), 
       ]

    f.writelines(L)
    f.close()

    arrays_to_save = [SPs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPRECE1,STELLARTYPEPRECE2,STELLARTYPEPSTCE1,STELLARTYPEPSTCE2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSPRECE1,MASSPRECE2,MASSPSTCE1,MASSPSTCE2,SEMIMAJORAXISSP,SEMIMAJORAXISPRECE, SEMIMAJORAXISPSTCE,MASKSPBH1,MASKSPBH2, 
                    MASKPRECEBH1,MASKPRECEBH2,MASKSPunb,MASKSPdco,MASKSPmrgr,MASKPRECEinSP,MASKPRECEBHNS1,MASKPRECEBHNS2, MASKPRECEnonBH, MASKPRECESEMAJ,
                    MASKPRECEnegCE, MASKPRECEposCE,MASKFIRSTCE,
                    MASKPSTCEBH1,MASKPSTCEBH2,MASKSPunb,MASKSPdco,MASKSPmrgr,MASKPSTCEinSP,MASKPSTCEBHNS1,MASKPSTCEBHNS2, MASKPSTCEnonBH, MASKPSTCESEMAJ,
                    MASKPSTCEnegCE, MASKPSTCEposCE,MASKLASTCE,
                    ORBITALPERIODSP, ORBITALPERIODPRECE,ORBITALPERIODPSTCE,TIMECE, MASKPRECERBPER, COSIPRECE,MASKPSTCERBPER, COSIPSTCE,
                    MASKPRECESEARCHABILITY, MASKPSTCESEARCHABILITY]
    filenames = ["SPs", "CEs", "STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2", "STELLARTYPEPRECE1", "STELLARTYPEPRECE2", "STELLARTYPEPSTCE1", "STELLARTYPEPSTCE2",
                    "MASSZAMSSP1", "MASSZAMSSP2", "MASSPRECE1", "MASSPRECE2", "MASSPSTCE1", "MASSPSTCE2", "SEMIMAJORAXISSP", "SEMIMAJORAXISPRECE", "SEMIMAJORAXISPSTCE", "MASKSPBH1", "MASKSPBH2", 
                    "MASKPRECEBH1", "MASKPRECEBH2", "MASKSPunb", "MASKSPdco", "MASKSPmrgr", "MASKPRECEinSP", "MASKPRECEBHNS1", "MASKPRECEBHNS2", "MASKPRECEnonBH", "MASKPRECESEMAJ", 
                    "MASKPRECEnegCE", "MASKPRECEposCE", "MASKFIRSTCE",
                    "MASKPSTCEBH1", "MASKPSTCEBH2", "MASKSPunb", "MASKSPdco", "MASKSPmrgr", "MASKPSTCEinSP", "MASKPSTCEBHNS1", "MASKPSTCEBHNS2", 
                    "MASKPSTCEnonBH", "MASKPSTCESEMAJ", "MASKPSTCEnegCE", "MASKPSTCEposCE", "MASKLASTCE",
                    "ORBITALPERIODSP", "ORBITALPERIODPRECE", "ORBITALPERIODPSTCE", "TIMECE", "MASKPRECERBPER", "COSIPRECE", "MASKPSTCERBPER", "COSIPSTCE",
                    "MASKPRECESEARCHABILITY", "MASKPSTCESEARCHABILITY"]

    i = 0
    if not os.path.exists("/data/a.saricaoglu/Files/" + IMF + "_CE" +"/" +  str(s.strftime("%m.%d"))): 
        os.makedirs("/data/a.saricaoglu/Files/" + IMF + "_CE" +"/" +  str(s.strftime("%m.%d"))) 
    for f in arrays_to_save:
        name = filenames[i]
        np.savetxt("/data/a.saricaoglu/Files/" + IMF + "_CE" +"/" +  str(s.strftime("%m.%d")) + "/" + name + "_" + str(s.strftime("%m.%d")) + ".txt", f)
        i+=1

e = datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time,
      " \n Seconds per sys. ", duration/len(SPs))
