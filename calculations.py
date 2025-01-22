import numpy as np
from datetime import datetime

# def get_current_parameters():

#     SPs = []
#     MTs = []
#     CEs = []

#     EVENTSMT = []
#     EVOLUTIONSTATSP = []

#     STELLARTYPEZAMSSP1   =  []
#     STELLARTYPEZAMSSP2   =  []
#     STELLARTYPESP1   =  []
#     STELLARTYPESP2   =  []       
#     STELLARTYPEPREMT1   =  []
#     STELLARTYPEPREMT2   =  [] 
#     STELLARTYPEPSTMT1   =  []
#     STELLARTYPEPSTMT2   =  [] 

#     MASSZAMSSP1 = []
#     MASSZAMSSP2 = []
#     MASSPREMT1 = []
#     MASSPREMT2 = []
#     MASSPSTMT1 = []
#     MASSPSTMT2 = []

#     SEMIMAJORAXISZAMSSP = []
#     SEMIMAJORAXISPREMT = []
#     SEMIMAJORAXISPSTMT = []

#     ORBITALPERIODSP = []
#     ORBITALPERIODPREMT = []
#     ORBITALPERIODPSTMT= []

#     RADIUSPREMT1 = []
#     RADIUSPREMT2 = []
#     RADIUSPSTMT1 = []
#     RADIUSPSTMT2 = []

#     TIMEPREMT = []
#     TIMEPSTMT = []

#     CEAFTERMT = []

#     arrays = [SPs,MTs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
#                     MASSZAMSSP1,MASSZAMSSP2,MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
#                     TIMEPREMT, TIMEPSTMT, CEAFTERMT, EVOLUTIONSTATSP, EVENTSMT,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2]


#     return arrays

def orbital_period(m1, m2, semimajax):
    G =  39.4769264 # gravitational constant in AU^3 / (year^2 x Msun) 
    M = (np.asarray(m1) + np.asarray(m2))
    A = np.asarray(semimajax)
    orbitalPeriod = (np.sqrt((4 * np.pi**2 * A**3)) / np.sqrt(G * M)) * 365.25

    return orbitalPeriod

def orbital_inclination(r1,r2,semimajax):

    # COMPAS output gives R in units of Rsun, hence converting AU here.
    # Also we consider the larger r_i for the calculations since if there is a BH it will always have smaller r and we want the star radius.
    R = np.maximum(np.asarray(r1),np.asarray(r2)) * 0.00465047
    print(R)
    cosi = R/np.asarray(semimajax)

    return cosi

def searchability(orbital_inclination):

    inc = np.random.uniform(-1,1,len(orbital_inclination))
    searchability_index = [np.absolute(inc) <= orbital_inclination][0]

    return searchability_index

def envelope_ejection(r1,r2,postcesemimajax):
    R = np.asarray(r1) + np.asarray(r2)
    ejection = [R < np.asarray(postcesemimajax)][0]

    return ejection

def type_change_MT(preMT1, pstMT1, preMT2, pstMT2, start_time):
    maskchangeMT1 = np.zeros(len(preMT1), dtype='bool')
    maskchangeMT2 = np.zeros(len(preMT2), dtype='bool')
    s = datetime.now()

    f = open("/data/a.saricaoglu/Files/Kroupa_MT/" + "MT_typechangeLog1_" + str(s.strftime("%m.%d")) +  ".txt", "a")
    f.writelines(["\n","\n Run :", start_time])

    for i in range(0, len(preMT1)):
        if (preMT1[i] == pstMT1[i]):
            maskchangeMT1[i] = False
            L = ['\n StellarType(1)<MT : ', str(preMT1[i]),
                 '\n StellarType(1)>MT : ', str(pstMT1[i]),
                 '\n No change in stellar type during mass transfer event, maskchangeMT1[',str(i), '] = ', str(maskchangeMT1[i])]

        else:
            maskchangeMT1[i] = True
            L = ['\n StellarType(1)<MT : ', str(preMT1[i]),
                 '\n StellarType(1)>MT : ', str(pstMT1[i]),
                 '\n Change in stellar type during mass transfer event, maskchangeMT1[',str(i), '] = ', str(maskchangeMT1[i])]
    f.writelines(L)
    f.close()

    f = open("/data/a.saricaoglu/Files/Kroupa_MT/" + "MT_typechangeLog2_" + str(s.strftime("%m.%d")) +  ".txt", "a")
    f.writelines(["\n","\n Run :", start_time])

    for i in range(0, len(preMT2)):
        if (preMT2[i] == pstMT2[i]):
            maskchangeMT2[i] = False
            L = ['\n StellarType(2)<MT : ', str(preMT2[i]),
                 '\n StellarType(2)>MT : ', str(pstMT2[i]),
                 '\n No change in stellar type during mass transfer event, maskchangeMT2[',str(i), '] = ', str(maskchangeMT2[i])]

        else:
            maskchangeMT2[i] = True
            L = ['\n StellarType(2)<MT : ', str(preMT2[i]),
                 '\n StellarType(2)>MT : ', str(pstMT2[i]),
                 '\n Change in stellar type during mass transfer event, maskchangeMT2[',str(i), '] = ', str(maskchangeMT2[i])]
    f.writelines(L)
    f.close()

    return maskchangeMT1, maskchangeMT2

def type_change_CE(preCE1, pstCE1, preCE2, pstCE2, start_time):
    maskchangeCE1 = np.zeros(len(preCE1), dtype='bool')
    maskchangeCE2 = np.zeros(len(preCE2), dtype='bool')

    f = open("/data/a.saricaoglu/Files/Kroupa_CE/" + "CE_typechangeLog1_" + str(s.strftime("%m.%d")) +  ".txt", "a")
    f.writelines(["\n","\n Run :", start_time])

    for i in range(0, len(preCE1)):
        if (preCE1[i] == pstCE1[i]):
            maskchangeCE1[i] = False
            L = ['\n StellarType(1)<CE : ', str(preCE1[i]),
                 '\n StellarType(1)>CE : ', str(pstCE1[i]),
                 '\n No change in stellar type during common envelope event, maskchangeCE1[',str(i), '] = ', str(maskchangeCE1[i])]

        else:
            maskchangeCE1[i] = True
            L = ['\n StellarType(1)<CE : ', str(preCE1[i]),
                 '\n StellarType(1)>CE : ', str(pstCE1[i]),
                 '\n Change in stellar type during common envelope event, maskchangeCE1[',str(i), '] = ', str(maskchangeCE1[i])]
    f.writelines(L)
    f.close()

    f = open("/data/a.saricaoglu/Files/Kroupa_CE/" + "CE_typechangeLog2_" + str(s.strftime("%m.%d")) +  ".txt", "a")
    f.writelines(["\n","\n Run :", start_time])

    for i in range(0, len(preCE2)):
        if (preCE2[i] == pstCE2[i]):
            maskchangeCE2[i] = False
            L = ['\n StellarType(2)<CE : ', str(preCE2[i]),
                 '\n StellarType(2)>CE : ', str(pstCE2[i]),
                 '\n No change in stellar type during common envelope event, maskchangeCE2[',str(i), '] = ', str(maskchangeCE2[i])]

        else:
            maskchangeCE2[i] = True
            L = ['\n StellarType(2)<CE : ', str(preCE2[i]),
                 '\n StellarType(2)>CE : ', str(pstCE2[i]),
                 '\n Change in stellar type during common envelope event, maskchangeCE2[',str(i), '] = ', str(maskchangeCE2[i])]
    f.writelines(L)
    f.close()

    return maskchangeCE1, maskchangeCE2


def find_last_mt(MTs):

    masklastMT = np.zeros(len(MTs), dtype='bool')
    last = 0
    for i in range(0, len(MTs)):
        if (MTs[i] == last):
            masklastMT[i-1] = False
            masklastMT[i] = True
        else:
            masklastMT[i] = True
            last = MTs[i]

    return masklastMT

s = searchability(np.random.uniform(-1,1,10))

print(np.shape(s))
print(s)