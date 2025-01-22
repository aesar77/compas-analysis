#### PARENT SCRIPT: secundus_processng_081424.py
#----> this version captures the system just before and after mass transfer (MT) event.
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

    SEMIMAJORAXISSP = []
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

    

    # Boolen arrays for masking.
    MASKSPBH1 = np.array([], dtype='bool')
    MASKSPBH2 = np.array([], dtype='bool')

    MASKPREMT = np.array([], dtype='bool')
    MASKPSTMT = np.array([], dtype='bool')

    MASKSPunb = np.array([], dtype='bool')
    MASKSPdco = np.array([], dtype='bool')
    MASKSPmrgr = np.array([], dtype='bool')

    MASKPREMTBH1 = np.array([], dtype='bool')
    MASKPSTMTBH1 = np.array([], dtype='bool')

    MASKPREMTBH2 = np.array([], dtype='bool')
    MASKPSTMTBH2 = np.array([], dtype='bool')

    MASKPREMTinSP = np.array([], dtype='bool')
    MASKPSTMTinSP = np.array([], dtype='bool')

    MASKPREMTBHNS1 = np.array([], dtype='bool')
    MASKPSTMTBHNS1 = np.array([], dtype='bool')

    MASKPREMTBHNS2 = np.array([], dtype='bool')
    MASKPSTMTBHNS2 = np.array([], dtype='bool')

    MASKPREMTnonBH = np.array([], dtype='bool')
    MASKPSTMTnonBH = np.array([], dtype='bool')

    MASKPREMTSEMAJ = np.array([], dtype='bool')
    MASKPSTMTSEMAJ = np.array([], dtype='bool')

    MASKPREMTRBPER = np.array([], dtype='bool')
    MASKPSTMTRBPER = np.array([], dtype='bool')

    MASKPREMTnegCE = np.array([], dtype='bool')
    MASKPSTMTnegCE = np.array([], dtype='bool')

    MASKPREMTposCE = np.array([], dtype='bool')
    MASKPSTMTposCE = np.array([], dtype='bool')

    COSIPREMT = np.array([], dtype='bool')
    COSIPSTMT = np.array([], dtype='bool')

    MASKFIRSTMT = np.array([], dtype='bool')
    MASKLASTMT = np.array([], dtype='bool')

    MASKPREMTSEARCHABILITY = np.array([], dtype='bool')
    MASKPSTMTSEARCHABILITY = np.array([], dtype='bool')

    MASSTRANSFERHISTORY = np.array([], dtype='bool')

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

        semimajorAxisSP = SP['SemiMajorAxis@ZAMS'][()] 
        semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] 
        semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

        timepreMT = MT['Time<MT'][()]
        timepstMT = MT['Time>MT'][()]
        
        massTransferhistory = MT['CEE>MT'][()]

        radiuspreMT1 = MT['Radius(1)<MT'][()]
        radiuspstMT1 = MT['Radius(1)>MT'][()] 
        radiuspreMT2 = MT['Radius(2)<MT'][()]
        radiuspstMT2 = MT['Radius(2)>MT'][()]   

        # Chooses first MT events.
        maskfirstMT = (eventsMT == 1)
        # Chooses last MT events
        masklastMT = np.zeros(len(MT['SEED'][()]), dtype='bool')
        last = 0
        for i in range(0, len(MT['SEED'][()])):
            if (MT['SEED'][()][i] == last):
                masklastMT[i-1] = False
                masklastMT[i] = True
            else:
                masklastMT[i] = True
                last = MT['SEED'][()][i]


        orbitalPeriodSP = calc.orbital_period(massZamsSP1, massZamsSP2, semimajorAxisSP)
        orbitalPeriodpreMT = calc.orbital_period(masspreMT1, masspreMT2, semimajorAxispreMT)
        orbitalPeriodpstMT = calc.orbital_period(masspstMT1, masspstMT2, semimajorAxispstMT)
        print(np.shape(orbitalPeriodpreMT))

        cosipreMT = calc.orbital_inclination(radiuspreMT1, radiuspreMT2, semimajorAxispreMT)
        cosipstMT = calc.orbital_inclination(radiuspstMT1, radiuspstMT2, semimajorAxispstMT)
        print(np.shape(cosipreMT))
        maskpreMTsearchability = calc.searchability(cosipreMT)
        maskpstMTsearchability = calc.searchability(cosipstMT)

        print(maskpreMTsearchability)
        print(np.shape(maskpreMTsearchability))
        # i = np.random.uniform(0,180,len(seedsMT))



## may add eccentricity later on

        maskSPBH1 = ((stellarTypeSP1  == 14) & (stellarTypeSP2  != 14)) 
        maskSPBH2 = ((stellarTypeSP1  != 14) & (stellarTypeSP2  == 14))
        maskpreMTBH1 = ((stellarTypepreMT1  == 14) & (stellarTypepreMT2  != 14)) 
        maskpreMTBH2 = ((stellarTypepreMT1  != 14) & (stellarTypepreMT2  == 14))       
        maskpstMTBH1 = ((stellarTypepstMT1  == 14) & (stellarTypepstMT2  != 14)) 
        maskpstMTBH2 = ((stellarTypepstMT1  != 14) & (stellarTypepstMT2  == 14))

        maskSPunb = (statusSP == 14)
        maskSPdco = (statusSP == 11)
        maskSPmrgr = (statusSP == 9) 

        maskMTinSP =  np.in1d(seedsSP, seedsMT)

        maskpreMTnonBH = ((stellarTypepreMT1 != 14) & (stellarTypepreMT2 != 14)) #none is bh
        maskpreMTBHNS1 = ((stellarTypepreMT1 == 14) & (stellarTypepreMT2 == 13)) #primary is bh secondary is ns
        maskpreMTBHNS2 = ((stellarTypepreMT1 == 13) & (stellarTypepreMT2 == 14)) #primary is ns secondary is bh

        maskpstMTnonBH = ((stellarTypepstMT1 != 14) & (stellarTypepstMT2 != 14)) #none is bh
        maskpstMTBHNS1 = ((stellarTypepstMT1 == 14) & (stellarTypepstMT2 == 13)) #primary is bh secondary is ns
        maskpstMTBHNS2 = ((stellarTypepstMT1 == 13) & (stellarTypepstMT2 == 14)) #primary is ns secondary is bh

        maskpreMTsemaj = (semimajorAxispreMT <= 3.0)
        maskpstMTsemaj = (semimajorAxispstMT <= 3.0)

        maskpreMTorbper = (orbitalPeriodpreMT <= 30)
        maskpstMTorbper = (orbitalPeriodpstMT <= 30)
        
        maskpreMTnegCE =  np.in1d(seedsMT*maskfirstMT, seedsCE, invert=True)
        maskpreMTposCE =  np.in1d(seedsMT*maskfirstMT, seedsCE)
        maskpstMTnegCE =  np.in1d(seedsMT*masklastMT, seedsCE, invert=True)
        maskpstMTposCE =  np.in1d(seedsMT*masklastMT, seedsCE)

        SPs.extend(seedsSP)
        MTs.extend(seedsMT)

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

        SEMIMAJORAXISSP.extend(semimajorAxisSP)
        SEMIMAJORAXISPREMT.extend(semimajorAxispreMT)
        SEMIMAJORAXISPSTMT.extend(semimajorAxispstMT)

        ORBITALPERIODSP.extend(orbitalPeriodSP)
        ORBITALPERIODPREMT.extend(orbitalPeriodpreMT)
        ORBITALPERIODPSTMT.extend(orbitalPeriodpstMT)

        MASSTRANSFERHISTORY = np.concatenate([MASSTRANSFERHISTORY, massTransferhistory])

        MASKSPunb = np.concatenate([MASKSPunb, maskSPunb])
        MASKSPdco = np.concatenate([MASKSPdco, maskSPdco])
        MASKSPmrgr = np.concatenate([MASKSPmrgr, maskSPmrgr])
        MASKMTinSP = np.concatenate([MASKPREMTinSP, maskMTinSP])


        MASKSPBH1 = np.concatenate([MASKSPBH1, maskSPBH1])
        MASKSPBH2 = np.concatenate([MASKSPBH2, maskSPBH2])
        MASKPREMTBH1 =np.concatenate([MASKPREMTBH1, maskpreMTBH1])
        MASKPREMTBH2 = np.concatenate([MASKPREMTBH2, maskpreMTBH2])

        MASKPREMTBHNS1 = np.concatenate([MASKPREMTBHNS1, maskpreMTBHNS1])
        MASKPREMTBHNS2 = np.concatenate([MASKPREMTBHNS2, maskpreMTBHNS2])

        MASKPREMTnonBH = np.concatenate([MASKPREMTnonBH, maskpreMTnonBH])
        MASKPREMTSEMAJ =  np.concatenate([MASKPREMTSEMAJ, maskpreMTsemaj])
        MASKPREMTRBPER = np.concatenate([MASKPREMTRBPER, maskpreMTorbper])
        COSIPREMT = np.concatenate([COSIPREMT, cosipreMT])
        MASKPREMTnegCE =  np.concatenate([MASKPREMTnegCE, maskpreMTnegCE])
        MASKPREMTposCE = np.concatenate([MASKPREMTposCE, maskpreMTposCE])

        MASKPSTMTBH1 =np.concatenate([MASKPSTMTBH1, maskpstMTBH1])
        MASKPSTMTBH2 = np.concatenate([MASKPSTMTBH2, maskpstMTBH2])

        MASKPSTMTBHNS1 = np.concatenate([MASKPSTMTBHNS1, maskpstMTBHNS1])
        MASKPSTMTBHNS2 = np.concatenate([MASKPSTMTBHNS2, maskpstMTBHNS2])

        MASKPSTMTnonBH = np.concatenate([MASKPSTMTnonBH, maskpstMTnonBH])
        MASKPSTMTSEMAJ =  np.concatenate([MASKPSTMTSEMAJ, maskpstMTsemaj])
        MASKPSTMTRBPER = np.concatenate([MASKPSTMTRBPER, maskpstMTorbper])
        COSIPSTMT = np.concatenate([COSIPSTMT, cosipstMT])
        MASKPSTMTnegCE =  np.concatenate([MASKPSTMTnegCE, maskpstMTnegCE])
        MASKPSTMTposCE = np.concatenate([MASKPSTMTposCE, maskpstMTposCE])

        MASKFIRSTMT = np.concatenate([MASKFIRSTMT, maskfirstMT])
        MASKLASTMT = np.concatenate([MASKLASTMT, masklastMT])

        MASKPREMTSEARCHABILITY = np.concatenate([MASKPREMTSEARCHABILITY, maskpreMTsearchability])
        MASKPSTMTSEARCHABILITY = np.concatenate([MASKPSTMTSEARCHABILITY, maskpstMTsearchability])

        c = datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 360 000 sys.)" + str(i) +  " end time :", current_time)
        i = i+1
    Data.close()

    L = ["\n Number of systems: ",str(len(MASSZAMSSP1)),  "\n Number of systems forming MT (via SP status): ", str(np.sum(MASKSPdco)), "\n Number of systems forming MT (via MT mass): ", str(len(MASSPREMT1*MASKFIRSTMT)), 
        "\n Number of primary mass BHs (SP): ", str(np.sum(MASKSPBH1)), "\n Number of primary mass BHs (MT): ", str(np.sum(MASKPREMTBH1)),  
        "\n Number of secondary mass BHs (SP): ", str(np.sum(MASKSPBH2)), "\n Number of secondary mass BHs (MT): ", str(np.sum(MASKPREMTBH2)),     
        "\n Number of unbound binaries (SP): ", str(np.sum(MASKSPunb)), 
        "\n Number of double compact objects (SP): ", str(np.sum(MASKSPdco)),  
        "\n Number of mergers (SP): ", str(np.sum(MASKSPmrgr)),
        "\n Number of BH-NS binaries (pre MT): ", str(np.sum(MASKPREMTBHNS1)), "\n Number of NS-BH binaries (pre MT): ", str(np.sum(MASKPREMTBHNS2)),
        "\n Number of systems without any BHs (pre MT): ", str(np.sum(MASKPREMTnonBH)),
        "\n Number of systems with semimajor axis < 3.0 AU (pre MT): ", str(np.sum(MASKPREMTSEMAJ)),
        "\n Number of systems undergoing CE (pre MT): ", str(np.sum(MASKPREMTposCE)),
        "\n Number of systems not undergoing CE (pre MT): ", str(np.sum(MASKPREMTnegCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEMAJ*MASKPREMTposCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND NOT undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEMAJ*MASKPREMTnegCE)),
        "\n Number of systems with orbital period < 30 days AND undergoing CE (pre MT): ", str(np.sum(MASKPREMTRBPER*MASKPREMTposCE)),
        "\n Number of systems with orbital period < 30 days AND NOT undergoing CE (pre MT): ", str(np.sum(MASKPREMTRBPER*MASKPREMTnegCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEARCHABILITY*MASKPREMTposCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND NOT undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEARCHABILITY*MASKPREMTnegCE)),
        "\n Number of BH-NS binaries (pst MT): ", str(np.sum(MASKPSTMTBHNS1)), "\n Number of NS-BH binaries (pst MT): ", str(np.sum(MASKPSTMTBHNS2)),
        "\n Number of systems without any BHs (pst MT): ", str(np.sum(MASKPSTMTnonBH)),
        "\n Number of systems with semimajor axis < 3.0 AU (pst MT): ", str(np.sum(MASKPSTMTSEMAJ)),
        "\n Number of systems undergoing CE (pst MT): ", str(np.sum(MASKPSTMTposCE)),
        "\n Number of systems not undergoing CE (pst MT): ", str(np.sum(MASKPSTMTnegCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEMAJ*MASKPSTMTposCE)),
        "\n Number of systems with semimajor axis < 3.0 AU AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEMAJ*MASKPSTMTnegCE)),
        "\n Number of systems with orbital period < 30 days AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTRBPER*MASKPSTMTposCE)),
        "\n Number of systems with orbital period < 30 days AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTRBPER*MASKPSTMTnegCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTposCE)),
        "\n Number of systems with orbital inclination (cosi) <  2.228e-6 AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTnegCE)),
        "\n Total MT events: ", str(np.sum(eventsMT)),
        "\n Number of primary HeWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT1 == 10 ])),
        "\n Number of primary COWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT1 == 11 ])),
        "\n Number of primary ONeWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT1 == 12 ])),
        "\n Number of primary HeWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT1 == 10 ])),
        "\n Number of primary COWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT1 == 11 ])),
        "\n Number of primary ONeWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT1 == 12 ])),
        "\n Number of secondary HeWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT2 == 10 ])),
        "\n Number of secondary COWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT2 == 11 ])),
        "\n Number of secondary ONeWD (pre MT) : ", str(np.sum([STELLARTYPEPREMT2 == 12 ])),
        "\n Number of secondary HeWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT2 == 10 ])),
        "\n Number of secondary COWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT2 == 11 ])),
        "\n Number of secondary ONeWD (pst MT) : ", str(np.sum([STELLARTYPEPSTMT2 == 12 ])), 

        ]

    f.writelines(L)
    f.close()

    arrays_to_save = [SPs,MTs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISSP,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,MASKSPBH1,MASKSPBH2, 
                    MASKPREMTBH1,MASKPREMTBH2,MASKSPunb,MASKSPdco,MASKSPmrgr,MASKPREMTinSP,MASKPREMTBHNS1,MASKPREMTBHNS2, MASKPREMTnonBH, MASKPREMTSEMAJ,
                    MASKPREMTnegCE, MASKPREMTposCE,MASKFIRSTMT,
                    MASKPSTMTBH1,MASKPSTMTBH2,MASKSPunb,MASKSPdco,MASKSPmrgr,MASKPSTMTinSP,MASKPSTMTBHNS1,MASKPSTMTBHNS2, MASKPSTMTnonBH, MASKPSTMTSEMAJ,
                    MASKPSTMTnegCE, MASKPSTMTposCE,MASKLASTMT,
                    ORBITALPERIODSP, ORBITALPERIODPREMT,ORBITALPERIODPSTMT,TIMEPREMT, TIMEPSTMT, MASKPREMTRBPER, COSIPREMT,MASKPSTMTRBPER, COSIPSTMT,
                    MASKPREMTSEARCHABILITY, MASKPSTMTSEARCHABILITY, MASSTRANSFERHISTORY]
    filenames = ["SPs", "MTs", "STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2", "STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                    "MASSZAMSSP1", "MASSZAMSSP2", "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISSP", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT", "MASKSPBH1", "MASKSPBH2", 
                    "MASKPREMTBH1", "MASKPREMTBH2", "MASKSPunb", "MASKSPdco", "MASKSPmrgr", "MASKPREMTinSP", "MASKPREMTBHNS1", "MASKPREMTBHNS2", "MASKPREMTnonBH", "MASKPREMTSEMAJ", 
                    "MASKPREMTnegCE", "MASKPREMTposCE", "MASKFIRSTMT",
                    "MASKPSTMTBH1", "MASKPSTMTBH2", "MASKSPunb", "MASKSPdco", "MASKSPmrgr", "MASKPSTMTinSP", "MASKPSTMTBHNS1", "MASKPSTMTBHNS2", 
                    "MASKPSTMTnonBH", "MASKPSTMTSEMAJ", "MASKPSTMTnegCE", "MASKPSTMTposCE", "MASKLASTMT",
                    "ORBITALPERIODSP", "ORBITALPERIODPREMT", "ORBITALPERIODPSTMT", "TIMEPREMT", "TIMEPSTMT", "MASKPREMTRBPER", "COSIPREMT", "MASKPSTMTRBPER", "COSIPSTMT",
                    "MASKPREMTSEARCHABILITY", "MASKPSTMTSEARCHABILITY","MASSTRANSFERHISTORY"]

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