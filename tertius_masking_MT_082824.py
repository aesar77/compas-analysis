#### PARENT SCRIPT: secundus_processng_081624_MT.py
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
IMFs = ["Kroupa_MT"]
date = "/09.02/"

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)

from compasUtils import printCompasDetails, getEventHistory, getEventStrings

# Choose an output hdf5 file to work with
for IMF in IMFs:

    pathToData = '/data/a.saricaoglu/Files/' + IMF + date #change the date accordingly the date of the files created via Sec
    # Keeps corresponding numerical values. 

    SPs = np.array([], dtype='float')
    MTs = np.array([], dtype='float')
    CEs = np.array([], dtype='float')

    EVENTSMT = np.array([], dtype='float')
    EVOLUTIONSTATSP = np.array([], dtype='float')

    STELLARTYPEZAMSSP1   = np.array([], dtype='float')
    STELLARTYPEZAMSSP2   = np.array([], dtype='float')
    STELLARTYPESP1   = np.array([], dtype='float')
    STELLARTYPESP2   = np.array([], dtype='float')       
    STELLARTYPEPREMT1   = np.array([], dtype='float')
    STELLARTYPEPREMT2   = np.array([], dtype='float') 
    STELLARTYPEPSTMT1   = np.array([], dtype='float')
    STELLARTYPEPSTMT2   = np.array([], dtype='float') 

    MASSZAMSSP1 = np.array([], dtype='float')
    MASSZAMSSP2 = np.array([], dtype='float')
    MASSPREMT1 = np.array([], dtype='float')
    MASSPREMT2 = np.array([], dtype='float')
    MASSPSTMT1 = np.array([], dtype='float')
    MASSPSTMT2 = np.array([], dtype='float')

    SEMIMAJORAXISZAMSSP = np.array([], dtype='float')
    SEMIMAJORAXISPREMT = np.array([], dtype='float')
    SEMIMAJORAXISPSTMT = np.array([], dtype='float')

    ORBITALPERIODSP = np.array([], dtype='float')
    ORBITALPERIODPREMT = np.array([], dtype='float')
    ORBITALPERIODPSTMT= np.array([], dtype='float')

    RADIUSPREMT1 = np.array([], dtype='float')
    RADIUSPREMT2 = np.array([], dtype='float')
    RADIUSPSTMT1 = np.array([], dtype='float')
    RADIUSPSTMT2 = np.array([], dtype='float')

    TIMEPREMT = np.array([], dtype='float')
    TIMEPSTMT = np.array([], dtype='float')

    

    # Boolen arrays for MASKing.
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

    CEAFTERMT = np.array([], dtype='bool')
    # Dataframes to use with sns plotting package
    data_SP = {}
    SPdf = pd.DataFrame(data_SP)
    data_MT = {}
    MTdf = pd.DataFrame(data_MT)
    data_CE = {}
    CEdf = pd.DataFrame(data_CE)
    
    # Reads files produced by ... and contructs the dataframes. SPs and MTs are separate since they have different sizes.
    arrays_to_save = [SPs,MTs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                    TIMEPREMT, TIMEPSTMT, CEAFTERMT, EVOLUTIONSTATSP, EVENTSMT,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2]
    filenames = ["SPs", "MTs","CEs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2", "STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                    "MASSZAMSSP1", "MASSZAMSSP2", "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISZAMSSP", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
                    "TIMEPREMT", "TIMEPSTMT","CEAFTERMT", "EVOLUTIONSTATSP", "EVENTSMT", "RADIUSPREMT1", "RADIUSPREMT2", "RADIUSPSTMT1", "RADIUSPSTMT2"]
    files = [x[2] for x in os.walk(pathToData)][0]
    i=0
    lenSP = 0
    lenMT = 0

    for file in files:
        print("file :", file)
        for var in filenames:
            print("var :", var)
            print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
            if var in file and len(arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]]) == 0 :
                index = np.where(np.asarray(filenames) == var)[0][0]
    

        print(pathToData  + file)
        data = np.loadtxt(pathToData  + file)
        print("index: ", index)
        print("File name: ", file, " Variable name: ", filenames[index])
        if len(data) != 0:
            print("data reached.")
        print("len: ", len(data))
        print("df col name:", filenames[index])   
        if type(arrays_to_save[index]) is list:
            arrays_to_save[index].append(data)
            print("type: list ---> SUCCESS")
            i = i+1
        else:
            arrays_to_save[index] = np.concatenate([arrays_to_save[index], data])
            print("type: ndarray ---> SUCCESS")

            i = i+1

    # for file in files:
    #     print("file :", file)
    #     for var in filenames:
    #         print("var :", var)
    #         print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
    #         if var in file and len(arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]]) == 0 :
    #             index = np.where(np.asarray(filenames) == var)[0][0]
    #             arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] = 1

    #     print(pathToData  + file)
    #     data = np.loadtxt(pathToData  + file)
    #     print("index: ", index)
    #     print("File name: ", file, " Variable name: ", filenames[index])
    #     if len(data) != 0:
    #         print("data reached.")
    #     print("len: ", len(data))
    #     print("df col name:", filenames[index])   

    #     if index == 0:
    #         lenSP = len(data)
    #     if index == 1:
    #         lenMT = len(data)
    #     if index == 2:
    #         lenCE = len(data_CE)

    #     if len(data) == lenSP:

    #         SPdf[filenames[index]] = data
    #         print(SPdf[filenames[index]] .isnull().sum(), SPdf[filenames[index]] .max())
    #         print("SUCCESS (SP)")
    #         i = i+1

    #     if len(data) == lenMT:

    #         MTdf[filenames[index]] = data
    #         print(MTdf[filenames[index]].isnull().sum(), MTdf[filenames[index]].max())
    #         print("SUCCESS (MT)")
    #         i = i+1

    #     if len(data) == lenCE:

    #         CEdf[filenames[index]] = data
    #         print(CEdf[filenames[index]].isnull().sum(), CEdf[filenames[index]].max())
    #         print("SUCCESS (CE)")
    #         i = i+1

    print(SPdf.keys())
    print(MTdf.keys())
    print(CEdf.keys())

    # Chooses first MT events.
    MASKFIRSTMT = (EVENTSMT == 1)

    # Chooses last MT events.
    MASKLASTMT = calc.find_last_mt(MTs)

    # Chooses MT events where stellar type changes.
    MASKTYPECHANGE1, MASKTYPECHANGE2 = calc.type_change_MT(STELLARTYPEPREMT1,STELLARTYPEPSTMT1,STELLARTYPEPREMT2, STELLARTYPEPSTMT2, start_time)


    ORBITALPERIODSP = calc.orbital_period(MASSZAMSSP1, MASSZAMSSP2, SEMIMAJORAXISZAMSSP)
    ORBITALPERIODPREMT = calc.orbital_period(MASSPREMT1, MASSPREMT2, SEMIMAJORAXISPREMT)
    ORBITALPERIODPSTMT = calc.orbital_period(MASSPSTMT2, MASSPSTMT2, SEMIMAJORAXISPSTMT)
    print(np.shape(ORBITALPERIODPREMT))

    COSIPREMT = calc.orbital_inclination(RADIUSPREMT1, RADIUSPREMT2, SEMIMAJORAXISPREMT)
    COSIPSTMT = calc.orbital_inclination(RADIUSPSTMT1, RADIUSPSTMT2, SEMIMAJORAXISPSTMT)
    print(np.shape(COSIPREMT))
    MASKPREMTSEARCHABILITY = calc.searchability(COSIPREMT)
    MASKPSTMTSEARCHABILITY = calc.searchability(COSIPSTMT)

    print(MASKPREMTSEARCHABILITY)
    print(np.shape(MASKPREMTSEARCHABILITY))
    # i = np.random.uniform(0,180,len(MTs))



## may add eccentricity later on

    MASKSPBH1 = ((STELLARTYPESP1  == 14) & (STELLARTYPESP2  != 14)) 
    MASKSPBH2 = ((STELLARTYPESP1  != 14) & (STELLARTYPESP2  == 14))
    MASKPREMTBH1 = ((STELLARTYPEPREMT1  == 14) & (STELLARTYPEPREMT2  != 14)) 
    MASKPREMTBH2 = ((STELLARTYPEPREMT1  != 14) & (STELLARTYPEPREMT2  == 14))       
    MASKPSTMTBH1 = ((STELLARTYPEPSTMT1  == 14) & (STELLARTYPEPSTMT2  != 14)) 
    MASKPSTMTBH2 = ((STELLARTYPEPSTMT1  != 14) & (STELLARTYPEPSTMT2  == 14))

    MASKSPunb = (EVOLUTIONSTATSP == 14)
    MASKSPdco = (EVOLUTIONSTATSP == 11)
    MASKSPmrgr = (EVOLUTIONSTATSP == 9) 

    MASKMTinSP =  np.in1d(SPs, MTs)

    MASKPREMTnonBH = ((STELLARTYPEPREMT1 != 14) & (STELLARTYPEPREMT2 != 14)) #none is bh
    MASKPREMTBHNS1 = ((STELLARTYPEPREMT1 == 14) & (STELLARTYPEPREMT2 == 13)) #primary is bh secondary is ns
    MASKPREMTBHNS2 = ((STELLARTYPEPREMT1 == 13) & (STELLARTYPEPREMT2 == 14)) #primary is ns secondary is bh

    MASKPSTMTnonBH = ((STELLARTYPEPSTMT1 != 14) & (STELLARTYPEPSTMT2 != 14)) #none is bh
    MASKPSTMTBHNS1 = ((STELLARTYPEPSTMT1 == 14) & (STELLARTYPEPSTMT2 == 13)) #primary is bh secondary is ns
    MASKPSTMTBHNS2 = ((STELLARTYPEPSTMT1 == 13) & (STELLARTYPEPSTMT2 == 14)) #primary is ns secondary is bh

    MASKPREMTsemaj = (SEMIMAJORAXISPREMT <= 3.0)
    MASKPSTMTsemaj = (SEMIMAJORAXISPSTMT <= 3.0)

    MASKPREMTorbper = (ORBITALPERIODPREMT <= 30)
    MASKPSTMTorbper = (ORBITALPERIODPSTMT <= 30)
    
    MASKPREMTnegCE =  np.in1d(MTs*MASKFIRSTMT, CEs, invert=True)
    MASKPREMTposCE =  np.in1d(MTs*MASKFIRSTMT, CEs)
    MASKPSTMTnegCE =  np.in1d(MTs*MASKLASTMT, CEs, invert=True)
    MASKPSTMTposCE =  np.in1d(MTs*MASKLASTMT, CEs)


f = open("/data/a.saricaoglu/Files/Kroupa/" + IMF + "_MT_" + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

f.writelines(["\n","\n Run :", start_time])

L = ["\n Number of systems: ",str(len(MASSZAMSSP1)),  "\n Number of systems forming MT (via SP STATUS): ", str(np.sum(MASKSPdco)), "\n Number of systems forming MT (via MT mass): ", str(len(MASSPREMT1*MASKFIRSTMT)), 
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
    "\n Number of systems with orbital inclination (COSI) <  2.228e-6 AND undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEARCHABILITY*MASKPREMTposCE)),
    "\n Number of systems with orbital inclination (COSI) <  2.228e-6 AND NOT undergoing CE (pre MT): ", str(np.sum(MASKPREMTSEARCHABILITY*MASKPREMTnegCE)),
    "\n Number of BH-NS binaries (pst MT): ", str(np.sum(MASKPSTMTBHNS1)), "\n Number of NS-BH binaries (pst MT): ", str(np.sum(MASKPSTMTBHNS2)),
    "\n Number of systems without any BHs (pst MT): ", str(np.sum(MASKPSTMTnonBH)),
    "\n Number of systems with semimajor axis < 3.0 AU (pst MT): ", str(np.sum(MASKPSTMTSEMAJ)),
    "\n Number of systems undergoing CE (pst MT): ", str(np.sum(MASKPSTMTposCE)),
    "\n Number of systems not undergoing CE (pst MT): ", str(np.sum(MASKPSTMTnegCE)),
    "\n Number of systems with semimajor axis < 3.0 AU AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEMAJ*MASKPSTMTposCE)),
    "\n Number of systems with semimajor axis < 3.0 AU AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEMAJ*MASKPSTMTnegCE)),
    "\n Number of systems with orbital period < 30 days AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTRBPER*MASKPSTMTposCE)),
    "\n Number of systems with orbital period < 30 days AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTRBPER*MASKPSTMTnegCE)),
    "\n Number of systems with orbital inclination (COSI) <  2.228e-6 AND undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTposCE)),
    "\n Number of systems with orbital inclination (COSI) <  2.228e-6 AND NOT undergoing CE (pst MT): ", str(np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTnegCE)),
    "\n Total MT events: ", str(np.sum(EVENTSMT)),
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

arrays_to_save = [MASKSPBH1,MASKSPBH2, MASKPREMTBH1,MASKPREMTBH2,MASKSPunb,MASKSPdco,MASKSPmrgr,MASKPREMTinSP,MASKPREMTBHNS1,MASKPREMTBHNS2, 
                MASKPREMTnonBH, MASKPREMTSEMAJ, MASKPREMTnegCE, MASKPREMTposCE,MASKFIRSTMT,
                MASKPSTMTBH1,MASKPSTMTBH2,MASKPSTMTinSP,MASKPSTMTBHNS1,MASKPSTMTBHNS2, MASKPSTMTnonBH, MASKPSTMTSEMAJ,
                MASKPSTMTnegCE, MASKPSTMTposCE, MASKLASTMT, MASKTYPECHANGE1, MASKTYPECHANGE2,
                ORBITALPERIODSP, ORBITALPERIODPREMT,ORBITALPERIODPSTMT, MASKPREMTRBPER, COSIPREMT,MASKPSTMTRBPER, COSIPSTMT,
                MASKPREMTSEARCHABILITY, MASKPSTMTSEARCHABILITY, CEAFTERMT]
filenames = [   "MASKSPBH1", "MASKSPBH2", "MASKPREMTBH1", "MASKPREMTBH2", "MASKSPunb", "MASKSPdco", "MASKSPmrgr", "MASKPREMTinSP", "MASKPREMTBHNS1", "MASKPREMTBHNS2", 
                "MASKPREMTnonBH", "MASKPREMTSEMAJ", "MASKPREMTnegCE", "MASKPREMTposCE", "MASKFIRSTMT",
                "MASKPSTMTBH1", "MASKPSTMTBH2", "MASKPSTMTinSP", "MASKPSTMTBHNS1", "MASKPSTMTBHNS2", "MASKPSTMTnonBH", "MASKPSTMTSEMAJ", 
                "MASKPSTMTnegCE", "MASKPSTMTposCE", "MASKLASTMT", "MASKTYPECHANGE1", "MASKTYPECHANGE2",
                "ORBITALPERIODSP", "ORBITALPERIODPREMT", "ORBITALPERIODPSTMT","MASKPREMTRBPER", "COSIPREMT", "MASKPSTMTRBPER", "COSIPSTMT",
                "MASKPREMTSEARCHABILITY", "MASKPSTMTSEARCHABILITY","CEAFTERMT"]

i = 0
if not os.path.exists("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d"))+ "/Masks/"): 
    os.makedirs("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d"))+ "/Masks/") 
for f in arrays_to_save:
    name = filenames[i]
    np.savetxt("/data/a.saricaoglu/Files/" + IMF + "_MT" +"/" +  str(s.strftime("%m.%d")) + "/Masks/" + name + "_" + str(s.strftime("%m.%d")) + ".txt", f)
    i+=1

e = datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time,
      " \n Seconds per sys. ", duration/len(SPs))