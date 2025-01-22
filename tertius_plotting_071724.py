####### PARENT SCRIPT: tertius_plotting_071724.py
#------> this version plots for contrained semimajor range < 0.5 AU for DCs
import os, sys
import numpy as np               # for handling arrays
import seaborn as sns
import numpy.ma as ma
import pandas as pd
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib
import seaborn as sns

sns.set_theme(style="ticks")

c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)

models = ["Kroupa"]
mode = "SingleBH"

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
# Choose an output hdf5 file to work with
for model in models:
    matplotlib.rcParams['figure.figsize'] = (21,14)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 17
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = '/data/a.saricaoglu/Files/' + model  + "/07.17/" #change the date accordingly the date of the files created via Sec

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

    # Boolean values for masking
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
   
    # Dataframes to use with sns plotting package
    data_SP = {}
    SPdf = pd.DataFrame(data_SP)
    data_DC = {}
    DCdf = pd.DataFrame(data_DC)
    
    # Reads files produced by ... and contructs the dataframes. SPs and DCs are separate since they have different sizes.
    arrays_to_save = [SPs,DCs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEDC1,STELLARTYPEDC2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSDC1,MASSDC2,SEMIMAJORAXISSP,SEMIMAJORAXISDC,MASKSPBH1,MASKSPBH2,MASKDCBH1,MASKDCBH2,
                    MASKSPunb,MASKSPdco,MASKSPmrgr,MASKDCinSP,MASKDCBHNS1,MASKDCBHNS2, MASKDCnonBH, MASKDCSEMAJ]
    filenames = ["SPs","DCs","STELLARTYPEZAMSSP1","STELLARTYPEZAMSSP2","STELLARTYPESP1","STELLARTYPESP2","STELLARTYPEDC1","STELLARTYPEDC2",
                "MASSZAMSSP1","MASSZAMSSP2","MASSDC1","MASSDC2","SEMIMAJORAXISSP","SEMIMAJORAXISDC","MASKSPBH1","MASKSPBH2","MASKDCBH1","MASKDCBH2",
                "MASKSPunb","MASKSPdco","MASKSPmrgr","MASKDCinSP","MASKDCBHNS1","MASKDCBHNS2", "MASKDCnonBH", "MASKDCSEMAJ"]

    runs= [x[2] for x in os.walk(pathToData)][0]
    i=0
    lenSP = 0
    lenDC = 0
    for run in runs:
        print("run :", run)
        for var in filenames:
            print("var :", var)
            print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
            if var in run and len(arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]]) == 0 :
                index = np.where(np.asarray(filenames) == var)[0][0]
                arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] = 1

        print(pathToData  + run)
        data = np.loadtxt(pathToData  + run)
        print("index: ", index)
        print("File name: ", run, " Variable name: ", filenames[index])
        if len(data) != 0:
            print("data reached.")
        print("len: ", len(data))
        print("df col name:", filenames[index])     
        if index == 0:
            lenSP = len(data)
        if index == 1:
            lenDC = len(data)

        if len(data) == lenSP:

            SPdf[filenames[index]] = data
            print(SPdf[filenames[index]] .isnull().sum(), SPdf[filenames[index]] .max())
            print("SUCCESS (SP)")
            i = i+1

        if len(data) == lenDC:

            DCdf[filenames[index]] = data
            print(DCdf[filenames[index]].isnull().sum(), DCdf[filenames[index]].max())
            print("SUCCESS (DC)")
            i = i+1


    print(SPdf.keys())
    print(DCdf.keys())

    systemSize = lenSP
    print("System size: ", lenSP, " 100/systemsize: ", 100/lenSP)
    DCdf_singleBH = pd.DataFrame()
    DCdf_tot = pd.DataFrame()
    #Systems which primary is BH
    DCdf_singleBH["BH1"] = DCdf["MASSDC1"]*DCdf["MASKDCBH1"]*DCdf["MASKDCSEMAJ"]
    DCdf_singleBH["CP2"] = DCdf["MASSDC2"]*DCdf["MASKDCBH1"]*DCdf["MASKDCSEMAJ"]
    DCdf_singleBH["SA1"] = DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCSEMAJ"]

    #Systems which secondary is BH
    DCdf_singleBH["BH2"] = DCdf["MASSDC2"]*DCdf["MASKDCBH2"]*DCdf["MASKDCSEMAJ"]
    DCdf_singleBH["CP1"] = DCdf["MASSDC1"]*DCdf["MASKDCBH2"]*DCdf["MASKDCSEMAJ"]
    DCdf_singleBH["SA2"] = DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCSEMAJ"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    DCdf_tot["BlackHole"] = np.concatenate([DCdf_singleBH["BH1"], DCdf_singleBH["BH2"]])
    DCdf_tot["Companion"] = np.concatenate([DCdf_singleBH["CP1"], DCdf_singleBH["CP2"]])
    DCdf_tot["Semax"] = np.concatenate([DCdf_singleBH["SA1"] , DCdf_singleBH["SA2"] ])
##########

    SPdf_singleBH = pd.DataFrame()
    SPdf_tot = pd.DataFrame()
     #Systems which primary is BH
    SPdf_singleBH["BH1"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH1"]
    SPdf_singleBH["CP2"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH1"]
    SPdf_singleBH["SA1"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH1"]

    #Systems which secondary is BH
    SPdf_singleBH["BH2"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH2"]
    SPdf_singleBH["CP1"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH2"]
    SPdf_singleBH["SA2"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH2"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    SPdf_tot["BlackHole"] = np.concatenate([SPdf_singleBH["BH1"], SPdf_singleBH["BH2"]])
    SPdf_tot["Companion"] = np.concatenate([SPdf_singleBH["CP1"], SPdf_singleBH["CP2"]])
    SPdf_tot["Semax"] = np.concatenate([SPdf_singleBH["SA1"] , SPdf_singleBH["SA2"] ])
##########
    
     #Systems which primary is BH (DCs in SP)
    SPdf_singleBH["BH1_DC"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["CP2_DC"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["SA1_DC"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]

    #Systems which secondary is BH (DCs in SP)
    SPdf_singleBH["BH2_DC"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["CP1_DC"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["SA2_DC"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary (DCs in SP)
    SPdf_tot["BlackHole_DC"] = np.concatenate([SPdf_singleBH["BH1_DC"], SPdf_singleBH["BH2_DC"]])
    SPdf_tot["Companion_DC"] = np.concatenate([SPdf_singleBH["CP1_DC"], SPdf_singleBH["CP2_DC"]])
    SPdf_tot["Semax_DC"] = np.concatenate([SPdf_singleBH["SA1_DC"] , SPdf_singleBH["SA2_DC"] ])
##########

    print(max(DCdf_tot["Semax"]))
    print(max(DCdf_tot["BlackHole"]))
    print(min(DCdf_tot["Semax"]))
    print(min(DCdf_tot["BlackHole"]))
    
    def mylog(x,y):
        return np.log10(x[y > 0]), np.log10(y[y > 0])
    def percentage(size, vals):
        return vals * (100/size)
    
######## UPDATE 07.11.24
    

    if not os.path.exists("/data/a.saricaoglu/Plots/" +  str(c.strftime("%m.%d"))): 
        os.makedirs("/data/a.saricaoglu/Plots/"  +  str(c.strftime("%m.%d"))) 

    # Produces heatmaps. Only black holes, companions are commented since not needed.
    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"],SPdf_tot["BlackHole"]),bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_DC"],SPdf_tot["BlackHole_DC"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_DC": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax_DC": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax_DC", columns="BlackHole_DC", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs in SP)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Mp_DC_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"],SPdf_tot["Companion"]), bins=20)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",  fontsize=40, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_DC"],SPdf_tot["Companion_DC"]), bins=20)
    # df = pd.DataFrame({
    # "Companion_DC": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax_DC": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion_DC", columns="Semax_DC", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs in SP)",  fontsize=40, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Ms_DC_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["Semax"],DCdf_tot["BlackHole"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems(DCO)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()



    # values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["Semax_DC"],DCdf_tot["Companion"]), bins=20)
    # values = percentage(systemSize, values)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DC)",  fontsize=40, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()


    #####

    