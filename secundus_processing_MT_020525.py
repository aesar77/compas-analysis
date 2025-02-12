#### PARENT SCRIPT: secundus_processng_082824_MT.py
#----> this version, instead of the systems just before the first and last mass transfer (MT) event, marks MT events in which StellarType1 or StellarType2
#      and gets the system parameters before and after each event. CAUTION: since each transition event for the system is stored as another (with the same seed number)
#      seed, arrays will be larger than actual system number.
import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time as t                      # for finding computation time
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
    MTs = []
    CEs = []

    EVENTSMT = []
    EVOLUTIONSTATSP = []

    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []       
    STELLARTYPEPREMT1   =  []
    STELLARTYPEPREMT2   =  [] 
    STELLARTYPEPSTMT1   =  []
    STELLARTYPEPSTMT2   =  [] 

    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    MASSPREMT1 = []
    MASSPREMT2 = []
    MASSPSTMT1 = []
    MASSPSTMT2 = []

    SEMIMAJORAXISZAMSSP = []
    SEMIMAJORAXISPREMT = []
    SEMIMAJORAXISPSTMT = []

    ORBITALPERIODSP = []
    ORBITALPERIODPREMT = []
    ORBITALPERIODPSTMT= []

    RADIUSPREMT1 = []
    RADIUSPREMT2 = []
    RADIUSPSTMT1 = []
    RADIUSPSTMT2 = []

    TIMEPREMT = []
    TIMEPSTMT = []

    CEAFTERMT = []


    f = open("/data/a.saricaoglu/Files/Kroupa/" + IMF + "_MT_" + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    f.writelines(["\n","\n Run :", start_time])

    i = 0

    for Data in data_outputs:

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " start time :", current_time)

        SP = Data['BSE_System_Parameters']
        MT = Data['BSE_RLOF']
        CE = Data['BSE_Common_Envelopes']

        seedsSP = SP['SEED'][()]
        seedsMT = MT['SEED'][()]
        seedsCE = CE['SEED'][()]

        eventsMT = MT['MT_Event_Counter'][()]

        statusSP = SP['Evolution_Status'][()]

        stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
        stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]        
        stellarTypepreMT1   =  MT['Stellar_Type(1)<MT'][()]
        stellarTypepreMT2   =  MT['Stellar_Type(2)<MT'][()]  
        stellarTypepstMT1   =  MT['Stellar_Type(1)>MT'][()]
        stellarTypepstMT2   =  MT['Stellar_Type(2)>MT'][()] 
        

        massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
        massZamsSP2 = SP['Mass@ZAMS(2)'][()] 

        # PreMT refers to the mass just before the first MT event, while pstMT refers to the mass just after the last MT event.
        masspreMT1 = MT['Mass(1)<MT'][()] 
        masspreMT2 = MT['Mass(2)<MT'][()] 
        masspstMT1 = MT['Mass(1)>MT'][()] 
        masspstMT2 = MT['Mass(2)>MT'][()]

        semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] 
        semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] 
        semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

        timepreMT = MT['Time<MT'][()]
        timepstMT = MT['Time>MT'][()]
        
        CEafterMT = MT['CEE>MT'][()]

        radiuspreMT1 = MT['Radius(1)<MT'][()]
        radiuspstMT1 = MT['Radius(1)>MT'][()] 
        radiuspreMT2 = MT['Radius(2)<MT'][()]
        radiuspstMT2 = MT['Radius(2)>MT'][()]   

        SPs.extend(seedsSP)
        MTs.extend(seedsMT)
        CEs.extend(seedsCE)

        EVENTSMT.extend(eventsMT)
        EVOLUTIONSTATSP.extend(statusSP)

        STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
        STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
        STELLARTYPESP1.extend(stellarTypeSP1)
        STELLARTYPESP2.extend(stellarTypeSP2)       
        STELLARTYPEPREMT1.extend(stellarTypepreMT1)
        STELLARTYPEPREMT2.extend(stellarTypepreMT2) 
        STELLARTYPEPSTMT1.extend(stellarTypepstMT1)
        STELLARTYPEPSTMT2.extend(stellarTypepstMT2) 

        MASSZAMSSP1.extend(massZamsSP1)
        MASSZAMSSP2.extend(massZamsSP2)
        MASSPREMT1.extend(masspreMT1)
        MASSPREMT2.extend(masspreMT2)
        MASSPSTMT1.extend(masspstMT1)
        MASSPSTMT2.extend(masspstMT2)

        SEMIMAJORAXISZAMSSP.extend(semimajorAxisZamsSP)
        SEMIMAJORAXISPREMT.extend(semimajorAxispreMT)
        SEMIMAJORAXISPSTMT.extend(semimajorAxispstMT)

        CEAFTERMT = np.concatenate([CEAFTERMT, CEafterMT])

        RADIUSPREMT1.extend(radiuspreMT1)
        RADIUSPREMT2.extend(radiuspreMT2)
        RADIUSPSTMT1.extend(radiuspstMT1)
        RADIUSPSTMT2.extend(radiuspstMT2)

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " end time :", current_time)
        i = i+1
    Data.close()

    arrays_to_save = [SPs,MTs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                    TIMEPREMT, TIMEPSTMT, CEAFTERMT, EVOLUTIONSTATSP, EVENTSMT]
    filenames = ["SPs", "MTs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2", "STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                    "MASSZAMSSP1", "MASSZAMSSP2", "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISZAMSSP", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
                    "TIMEPREMT", "TIMEPSTMT","CEAFTERMT", "EVOLUTIONSTATSP", "EVENTSMT"]

    i = 0
    if not os.path.exists("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d"))): 
        os.makedirs("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d"))) 
    for f in arrays_to_save:
        name = filenames[i]
        np.savetxt("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d")) + "/" + name + "_" + str(s.strftime("%m.%d")) + ".txt", f)
        i+=1

e = datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time,
      " \n Seconds per sys. ", duration/len(SPs))