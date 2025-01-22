#### PARENT SCRIPT: secundus_processng_MT_081624_MT.py
#----> this version, instead of only MT data frame, reads CE too.
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
    # SPs
    SPs = []
    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []   
    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    EVOLUTIONSTATSP = []
    SEMIMAJORAXISZAMSSP = []


    MTs = []
    EVENTSMT = []
    STELLARTYPEPREMT1   =  []
    STELLARTYPEPREMT2   =  [] 
    STELLARTYPEPSTMT1   =  []
    STELLARTYPEPSTMT2   =  [] 
    MASSPREMT1 = []
    MASSPREMT2 = []
    MASSPSTMT1 = []
    MASSPSTMT2 = []
    SEMIMAJORAXISPREMT = []
    SEMIMAJORAXISPSTMT = []
    RADIUSPREMT1 = []
    RADIUSPREMT2 = []
    RADIUSPSTMT1 = []
    RADIUSPSTMT2 = []
    TIMEPREMT = []
    TIMEPSTMT = []

    CEs = []
    EVENTSCE = []
    STELLARTYPEPRECE1   =  []
    STELLARTYPEPRECE2   =  [] 
    STELLARTYPEPSTCE1   =  []
    STELLARTYPEPSTCE2   =  [] 
    MASSPRECE1 = []
    MASSPRECE2 = []
    MASSPSTCE1 = []
    MASSPSTCE2 = []
    SEMIMAJORAXISPRECE = []
    SEMIMAJORAXISPSTCE = []
    RADIUSPRECE1 = []
    RADIUSPRECE2 = []
    RADIUSPSTCE1 = []
    RADIUSPSTCE2 = []
    COSIPRECE = []
    COSIPSTCE = []
    TIMECE = []

    


    f = open("/data/a.saricaoglu/Files/Kroupa/" + IMF + "_MT_" + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    f.writelines(["\n","\n Run :", start_time])

    i = 0

    for Data in data_outputs:

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " start time :", current_time)

        SP = Data['BSE_System_Parameters']

        seedsSP = SP['SEED'][()]
        statusSP = SP['Evolution_Status'][()]

        stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
        stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]  

        massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
        massZamsSP2 = SP['Mass@ZAMS(2)'][()]

        semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] 

        SPs.extend(seedsSP)
        EVOLUTIONSTATSP.extend(statusSP)

        STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
        STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
        STELLARTYPESP1.extend(stellarTypeSP1)
        STELLARTYPESP2.extend(stellarTypeSP2)   

        MASSZAMSSP1.extend(massZamsSP1)
        MASSZAMSSP2.extend(massZamsSP2)

        SEMIMAJORAXISZAMSSP.extend(semimajorAxisZamsSP)

        MT = Data['BSE_RLOF']

        seedsMT = MT['SEED'][()]
        eventsMT = MT['MT_Event_Counter'][()]

        stellarTypepreMT1   =  MT['Stellar_Type(1)<MT'][()]
        stellarTypepreMT2   =  MT['Stellar_Type(2)<MT'][()]  
        stellarTypepstMT1   =  MT['Stellar_Type(1)>MT'][()]
        stellarTypepstMT2   =  MT['Stellar_Type(2)>MT'][()] 

        masspreMT1 = MT['Mass(1)<MT'][()] 
        masspreMT2 = MT['Mass(2)<MT'][()] 
        masspstMT1 = MT['Mass(1)>MT'][()] 
        masspstMT2 = MT['Mass(2)>MT'][()]

        semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] 
        semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

        timepreMT = MT['Time<MT'][()]
        timepstMT = MT['Time>MT'][()]
        
        massTransferhistory = MT['CEE>MT'][()]

        radiuspreMT1 = MT['Radius(1)<MT'][()]
        radiuspstMT1 = MT['Radius(1)>MT'][()] 
        radiuspreMT2 = MT['Radius(2)<MT'][()]
        radiuspstMT2 = MT['Radius(2)>MT'][()]

        MTs.extend(seedsMT)
        EVENTSMT.extend(eventsMT)
      
        STELLARTYPEPREMT1.extend(stellarTypepreMT1)
        STELLARTYPEPREMT2.extend(stellarTypepreMT2) 
        STELLARTYPEPSTMT1.extend(stellarTypepstMT1)
        STELLARTYPEPSTMT2.extend(stellarTypepstMT2) 

        MASSPREMT1.extend(masspreMT1)
        MASSPREMT2.extend(masspreMT2)
        MASSPSTMT1.extend(masspstMT1)
        MASSPSTMT2.extend(masspstMT2)

        SEMIMAJORAXISPREMT.extend(semimajorAxispreMT)
        SEMIMAJORAXISPSTMT.extend(semimajorAxispstMT)

        TIMEPREMT.extend(timepreMT)
        TIMEPSTMT.extend(timepstMT)

        MASSTRANSFERHISTORY = np.concatenate([massTransferhistory])

        CE = Data['BSE_Common_Envelopes']

        seedsCE = CE['SEED'][()]
        eventsCE = CE['CE_Event_Counter'][()]
       
        stellarTypepreCE1   =  CE['Stellar_Type(1)<CE'][()]
        stellarTypepreCE2   =  CE['Stellar_Type(2)<CE'][()]  
        stellarTypepstCE1   =  CE['Stellar_Type(1)'][()]
        stellarTypepstCE2   =  CE['Stellar_Type(2)'][()] 
        # PreCE refers to the mass just before the first CE event, while pstCE refers to the mass just after the last CE event.
        masspreCE1 = CE['Mass(1)<CE'][()] 
        masspreCE2 = CE['Mass(2)<CE'][()] 
        masspstCE1 = CE['Mass(1)>CE'][()] 
        masspstCE2 = CE['Mass(2)>CE'][()]

        semimajorAxispreCE = CE['SemiMajorAxis<CE'][()] 
        semimajorAxispstCE = CE['SemiMajorAxis>CE'][()]

        timeCE = CE['Time'][()]
        
        radiuspreCE1 = CE['Radius(1)<CE'][()]
        radiuspstCE1 = CE['Radius(1)>CE'][()] 
        radiuspreCE2 = CE['Radius(2)<CE'][()]
        radiuspstCE2 = CE['Radius(2)>CE'][()] 

        CEs.extend(seedsCE)
        EVENTSCE.extend(eventsCE)
   
        STELLARTYPEPRECE1.extend(stellarTypepreCE1)
        STELLARTYPEPRECE2.extend(stellarTypepreCE2) 
        STELLARTYPEPSTCE1.extend(stellarTypepstCE1)
        STELLARTYPEPSTCE2.extend(stellarTypepstCE2) 

        MASSPRECE1.extend(masspreCE1)
        MASSPRECE2.extend(masspreCE2)
        MASSPSTCE1.extend(masspstCE1)
        MASSPSTCE2.extend(masspstCE2)

        SEMIMAJORAXISPRECE.extend(semimajorAxispreCE)
        SEMIMAJORAXISPSTCE.extend(semimajorAxispstCE)

      

        
 

        # PreMT refers to the mass just before the first MT event, while pstMT refers to the mass just after the last MT event.
   

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

        MASSTRANSFERHISTORY = np.concatenate([MASSTRANSFERHISTORY, massTransferhistory])

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " end time :", current_time)
        i = i+1
    Data.close()

    arrays_to_save = [SPs,MTs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                    TIMEPREMT, TIMEPSTMT, MASSTRANSFERHISTORY, EVOLUTIONSTATSP, EVENTSMT]
    filenames = ["SPs", "MTs","CEs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2", "STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                    "MASSZAMSSP1", "MASSZAMSSP2", "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISZAMSSP", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
                    "TIMEPREMT", "TIMEPSTMT","MASSTRANSFERHISTORY", "EVOLUTIONSTATSP", "EVENTSMT"]

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