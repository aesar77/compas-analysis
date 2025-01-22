####### PARENT SCRIPT: seabornplots_2dhist_singleBH copy. py
#------> this version sums up primary and secondary single bhs as primary one is BH and secondary one is star.

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

#matplotlib.use("Agg")
c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)

models = ["Kroupa"]
mode = "SingleBH"
# plt.figure(figsize=(12, 12))

# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
# Choose an output hdf5 file to work with
for model in models:
    matplotlib.rcParams['figure.figsize'] = (15,10)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 17
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = '/data/a.saricaoglu/Files/' + model  + "/07.14/"


    DCs = []
    MTs = []
    CEs = []
    SNe = []
    DCs = []

    mP_innms = 0
    mP_ins = 0
    mP_preCEs = 0
    mP_pstCEs = 0
    mP_DCs = 0
    mP_DCs_CEs = 0

    mS_innms = 0
    mS_ins = 0
    mS_preCEs = 0
    mS_pstCEs = 0
    mS_DCs = 0
    mS_DCs_CEs =0

    semax_innms = 0
    semax_ins = 0
    semax_preCEs = 0
    semax_pstCEs = 0
    semax_DCs = 0
    semax_DCs_CEs = 0

    maskSP = 0
    # maskCEpre = np.array(0, dtype='bool')
    # maskCEpst = np.array(0, dtype='bool') #mask_unique
    maskCEMult = 0
    maskpreCEBH = 0
    maskpstCEBH = 0
    # maskpstCEBH =np.array(0, dtype='bool')
    maskpreCEBBH =0
    maskpstCEBBH =0

    maskSPBH = 0
    maskSPBHpr = 0
    maskSPBHsc = 0

    maskDCBH = 0
    maskDCBHpr = 0
    maskDCBHsc = 0

    maskSPBBH =0
    maskDCBBH =0


    maskSPDCO = 0
    maskSPUnbound =  0
    maskSPMerg = 0

    maskSPBH_CEpos = 0
    maskSPBHpr_CEpos = 0
    maskSPBHsc_CEpos = 0

    maskDCBH_CEpos = 0
    maskDCBHpr_CEpos = 0
    maskDCBHsc_CEpos = 0

    maskSPBBH_CEpos =0
    maskDCBBH_CEpos =0


    maskSPDCO_CEpos = 0
    maskSPUnbound_CEpos =  0
    maskSPMerg_CEpos = 0
    maskSPBoundBinary = 0
    maskSPBoundBinary_CEpos = 0

    data_SP = {}
    SPdf = pd.DataFrame(data_SP)
    data_CP = {}
    CPdf = pd.DataFrame(data_CP)
    data_DC = {}
    DCdf = pd.DataFrame(data_DC)
    
    arrays_to_save = [mP_innms, mP_ins, mP_preCEs, mP_pstCEs, mP_DCs, mP_DCs_CEs, mS_innms, mS_ins, mS_preCEs, mS_pstCEs, mS_DCs, mS_DCs_CEs, semax_innms, semax_ins, semax_preCEs, semax_pstCEs, semax_DCs, semax_DCs_CEs,
                      maskpreCEBH, maskpstCEBH, maskpreCEBBH, maskpstCEBBH, maskSPBH, maskSPBHpr, maskSPBHsc, maskSPBBH, maskDCBH, maskDCBHpr, maskDCBHsc, maskDCBBH, 
                      maskSPUnbound, maskSPDCO, maskSPMerg, 
                      maskSPBH_CEpos, maskSPBHpr_CEpos, maskSPBHsc_CEpos, maskSPBBH_CEpos, maskDCBH_CEpos, maskDCBHpr_CEpos, maskDCBHsc_CEpos, maskDCBBH_CEpos, maskSPBBH_CEpos, maskSPUnbound_CEpos, maskSPDCO_CEpos, maskSPMerg_CEpos, maskCEMult, maskSPBoundBinary, maskSPBoundBinary_CEpos]

    filenames =  ['mP_innms', 'mP_ins', 'mP_preCEs', 'mP_pstCEs', 'mP_DCs', 'mP_DCs_CEs', 'mS_innms', 'mS_ins', 'mS_preCEs', 'mS_pstCEs', 'mS_DCs', 'mS_DCs_CEs', 'semax_innms', 'semax_ins', 'semax_preCEs', 'semax_pstCEs', 'semax_DCs','semax_DCs_CEs', 
                  'maskpreCEBH', 'maskpstCEBH', 'maskpreCEBBH', 'maskpstCEBBH', 'maskSPBH', 'maskSPBHpr', 'maskSPBHsc', 'maskSPBBH', 'maskDCBH', 'maskDCBHpr', 'maskDCBHsc', 'maskDCBBH',
                  'maskSPUnbound', 'maskSPDCO', 'maskSPMerg', 
                  'maskSPBH_CEpos', 'maskSPBHpr_CEpos', 'maskSPBHsc_CEpos', 'maskSPBBH_CEpos', 'maskDCBH_CEpos', 'maskDCBHpr_CEpos', 'maskDCBHsc_CEpos', 'maskDCBBH_CEpos', 'maskSPBBH_CEpos', 'maskSPUnbound_CEpos', 'maskSPDCO_CEpos', 'maskSPMerg_CEpos', 'maskCEMult', 'maskSPBoundBinary', 'maskSPBoundBinary_CEpos']

    runs= [x[2] for x in os.walk(pathToData)][0]
    # print(runs)
    i=0
    lenSP = 0
    lenCP = 0
    lenDC = 0
    for run in runs:
        print("run :", run)
        #run = run.strip()
        for var in filenames:
            print("var :", var)
            print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
            if var in run and arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] == 0 :
                index = np.where(np.asarray(filenames) == var)[0][0]
                arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] = 1

        data = np.loadtxt(pathToData + "/" + run)
        print("index: ", index)
        print("File name: ", run, " Variable name: ", filenames[index])
        if len(data) != 0:
            print("data reached.")
        print("len: ", len(data))
        print("df col name:", filenames[index])     
        if index == 0:
            lenSP = len(data)
        if index == 2:
            lenCP = len(data)
        if index == 4:
            lenDC = len(data)

        if len(data) == lenSP:

            SPdf[filenames[index]] = data
            print(SPdf[filenames[index]] .isnull().sum(), SPdf[filenames[index]] .max())
            print("SUCCESS (SP)")
            # print(df[filenames[index]][0:10])
            i = i+1

        if len(data) == lenCP:

            CPdf[filenames[index]] = data
            print(CPdf[filenames[index]].isnull().sum(), CPdf[filenames[index]].max())
            print("SUCCESS (CP)")
            # print(df[filenames[index]][0:10])
            i = i+1

        if len(data) == lenDC:

            DCdf[filenames[index]] = data
            print(DCdf[filenames[index]].isnull().sum(), DCdf[filenames[index]].max())
            print("SUCCESS (DC)")
            # print(df[filenames[index]][0:10])
            i = i+1


    print(SPdf.keys())
    print(CPdf.keys())
    print(DCdf.keys())

    print("semaxpreces nonfloat check: ", np.isnan(CPdf["semax_preCEs"]), "min: ", min(CPdf["semax_preCEs"]), " max: ", max(CPdf["semax_preCEs"]))
    print("mpreces nonfloat check: ", np.isnan(CPdf["mP_preCEs"]), "min: ", min(CPdf["mP_preCEs"]), " max: ", max(CPdf["mP_preCEs"]))
    # g = sns.JointGrid(data=SPdf, x="semax_innms", y="mP_innms", marginal_ticks=True)
    # #g.ax_joint.set(yscale="log")
    # cax = g.figure.add_axes([.25, .8, .5, .015])
    # g.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #      binrange = ((np.log10(min(SPdf["semax_innms"])), np.log10(max(SPdf["semax_innms"]))),(np.log10(min(SPdf["mP_innms"])), np.log10(max(SPdf["mP_innms"])))), thresh=None)
    # g.plot_marginals(sns.histplot, element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # g.ax_joint.set(xlabel="Semimajor Axis [AU]")
    # g.ax_joint.set(ylabel="$M_p$  [$M_{\odot}$]")
    # g.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_1_" + current_time + ".png",   bbox_inches='tight')

    # h = sns.JointGrid(data=SPdf, x="semax_innms", y="mS_innms", marginal_ticks=True)
    # #h.ax_joint.set(yscale="log")
    # cax = h.figure.add_axes([.25, .8, .5, .015])
    # h.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #    binrange = (( np.log10(min(SPdf["semax_innms"])), np.log10(max(SPdf["semax_innms"]))),(np.log10(min(SPdf["mS_innms"])),np.log10(max(SPdf["mS_innms"])))), thresh=None)
    # h.plot_marginals(sns.histplot, element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mS_ins",levels=10, color="b", linewidths=1)
    # h.ax_joint.set(xlabel="Semimajor Axis [AU]")
    # h.ax_joint.set(ylabel="$M_s$  [$M_{\odot}$]")
    # h.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode + "_2d_hist_2_" + current_time + ".png",   bbox_inches='tight')  


    # print(((np.log10(0.1), np.log10(max(CPdf["semax_preCEs"]))),(np.log10(0.1),np.log10(max(CPdf["mP_preCEs"])))))
    # j = sns.JointGrid(data=CPdf, x="semax_preCEs", y="mP_preCEs", marginal_ticks=True)
    # #j.ax_joint.set(yscale="log")

    # cax = j.figure.add_axes([.25, .8, .5, .015])
    # j.plot_joint(
    #     sns.histplot, discrete=(False, False), log_scale = (True, True),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_preCEs"]))), (np.log10(0.1), np.log10(max(CPdf["mP_preCEs"])))), thresh=None)
    # j.plot_marginals(sns.histplot, binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_preCEs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)bin
    # j.ax_joint.set(xlabel="Semimajor Axis before CE [AU]")
    # j.ax_joint.set(ylabel="$M_p$ before CE [$M_{\odot}$]")
    # j.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode + "_2d_hist_3_" + current_time + ".png",   bbox_inches='tight')

    # k = sns.JointGrid(data=CPdf, x="semax_preCEs", y="mS_preCEs", marginal_ticks=True)
    # #k.ax_joint.set(yscale="log")
    # cax = k.figure.add_axes([.25, .8, .5, .015])
    # k.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #         binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_preCEs"]))), (np.log10(0.1), np.log10(max(CPdf["mS_preCEs"])))),  thresh=None)
    # k.plot_marginals(sns.histplot, binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_preCEs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # k.ax_joint.set(xlabel="Semimajor Axis before CE [AU]")
    # k.ax_joint.set(ylabel="$M_s$ before CE [$M_{\odot}$]")
    # k.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + model  + "_" +  mode + "_2d_hist_4_" + current_time + ".png",   bbox_inches='tight')


    # l = sns.JointGrid(data=CPdf, x="semax_pstCEs", y="mP_pstCEs", marginal_ticks=True)
    # #l.ax_joint.set(yscale="log")
    # cax = l.figure.add_axes([.25, .8, .5, .015])
    # l.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_pstCEs"]))), (np.log10(0.1), np.log10(max(CPdf["mP_pstCEs"])))),thresh=None)
    # l.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_pstCEs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # l.ax_joint.set(xlabel="Semimajor Axis  [AU]")
    # l.ax_joint.set(ylabel="$M_p$  [$M_{\odot}$]")
    # l.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + model  + "_" +  mode + "_2d_hist_5_" + current_time + ".png",   bbox_inches='tight')

    # m = sns.JointGrid(data=CPdf, x="semax_pstCEs", y="mS_pstCEs", marginal_ticks=True)
    # #m.ax_joint.set(yscale="log")
    # cax = m.figure.add_axes([.25, .8, .5, .015])
    # m.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_pstCEs"]))), (np.log10(0.1), np.log10(max(CPdf["mP_pstCEs"])))), thresh=None)
    # m.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(CPdf["semax_pstCEs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=SPdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # m.ax_joint.set(xlabel="Semimajor Axis  [AU]")
    # m.ax_joint.set(ylabel="$M_s$  [$M_{\odot}$]")
    # m.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_6_" + current_time + ".png",   bbox_inches='tight')



    # print(DCdf["mP_DCs"])
    # n = sns.JointGrid(data=DCdf, x="semax_DCs", y="mP_DCs", marginal_ticks=True)
    # #n.ax_joint.set(yscale="log")
    # cax = n.figure.add_axes([.25, .8, .5, .015])
    # n.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"]))), (np.log10(0.1), np.log10(max(DCdf["mP_DCs"])))), thresh=None)
    # n.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=DCdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # n.ax_joint.set(xlabel="Semimajor Axis DCO [AU]")
    # n.ax_joint.set(ylabel="$M_p$ DCO [$M_{\odot}$]")
    # n.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_7_" + current_time + ".png",   bbox_inches='tight')
    
    # print(DCdf["mS_DCs"])

    # n = sns.JointGrid(data=DCdf, x="semax_DCs", y="mS_DCs", marginal_ticks=True)
    # #n.ax_joint.set(yscale="log")
    # cax = n.figure.add_axes([.25, .8, .5, .015])
    # n.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"]))), (np.log10(0.1), np.log10(max(DCdf["mS_DCs"])))), thresh=None)
    # n.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=DCdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # n.ax_joint.set(xlabel="Semimajor Axis DCO [AU]")
    # n.ax_joint.set(ylabel="$M_s$ DCO [$M_{\odot}$]")
    # n.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_8_" + current_time + ".png",   bbox_inches='tight')

    # print(DCdf["mP_DCs_CEs"])

    # n = sns.JointGrid(data=DCdf, x="semax_DCs_CEs", y="mP_DCs_CEs", marginal_ticks=True)
    # #n.ax_joint.set(yscale="log")
    # cax = n.figure.add_axes([.25, .8, .5, .015])
    # n.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs_CEs"]))), (np.log10(0.1), np.log10(max(DCdf["mP_DCs_CEs"])))), thresh=None)
    # n.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=DCdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # n.ax_joint.set(xlabel="Semimajor Axis DCO CE positive [AU]")
    # n.ax_joint.set(ylabel="$M_p$ DCO CE positive [$M_{\odot}$]")
    # n.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_9_" + current_time + ".png",   bbox_inches='tight')

    # n = sns.JointGrid(data=DCdf, x="semax_DCs_CEs", y="mS_DCs_CEs", marginal_ticks=True)
    # #n.ax_joint.set(yscale="log")
    # cax = n.figure.add_axes([.25, .8, .5, .015])
    # n.plot_joint(
    #     sns.histplot, discrete=(False, False),
    #     cmap="BuPu", pmax=1.0, cbar=True, cbar_ax=cax, cbar_kws= {"orientation" :"horizontal"}, log_scale = (True, True), common_bins=False,
    #     binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs_CEs"]))), (np.log10(0.1), np.log10(max(DCdf["mS_DCs_CEs"])))), thresh=None)
    # n.plot_marginals(sns.histplot,binrange = ((np.log10(0.1), np.log10(max(DCdf["semax_DCs"])))), element="bars", color="#03012d")
    # # sns.kdeplot(data=DCdf, x="semax_ins", y="mP_ins",levels=10, color="b", linewidths=1)
    # n.ax_joint.set(xlabel="Semimajor Axis DCO CE positive [AU]")
    # n.ax_joint.set(ylabel="$M_s$ DCO CE positive[$M_{\odot}$]")
    # n.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d"))+ "/" + model  + "_" +  mode +  "_2d_hist_10_" + current_time + ".png",   bbox_inches='tight')


    DCdf_singleBH = pd.DataFrame()
    DCdf_tot = pd.DataFrame()
    DCdf_singleBH["mP"] = DCdf["mP_DCs"]*DCdf["maskDCBH"]
    DCdf_singleBH["mS"] = DCdf["mS_DCs"]*DCdf["maskDCBH"]
    DCdf_singleBH["semax"] = DCdf["semax_DCs"]*DCdf["maskDCBH"]

    DCdf_singleBH["mP_pr"] = DCdf["mP_DCs"]*DCdf["maskDCBHpr"]
    DCdf_singleBH["mS_pr"] = DCdf["mS_DCs"]*DCdf["maskDCBHpr"]
    DCdf_singleBH["semax_pr"] = DCdf["semax_DCs"]*DCdf["maskDCBHpr"]

    DCdf_singleBH["mP_sc"] = DCdf["mP_DCs"]*DCdf["maskDCBHsc"]
    DCdf_singleBH["mS_sc"] = DCdf["mS_DCs"]*DCdf["maskDCBHsc"]
    DCdf_singleBH["semax_sc"] = DCdf["semax_DCs"]*DCdf["maskDCBHsc"]

    DCdf_singleBH["mP_CE"] = DCdf["mP_DCs_CEs"]*DCdf["maskDCBH_CEpos"]
    DCdf_singleBH["mS_CE"] = DCdf["mS_DCs_CEs"]*DCdf["maskDCBH_CEpos"]
    DCdf_singleBH["semax_CE"] = DCdf["semax_DCs_CEs"]*DCdf["maskDCBH_CEpos"]

    DCdf_singleBH["mP_pr_CE"] = DCdf["mP_DCs_CEs"]*DCdf["maskDCBHpr_CEpos"]
    DCdf_singleBH["mS_pr_CE"] = DCdf["mS_DCs_CEs"]*DCdf["maskDCBHpr_CEpos"]
    DCdf_singleBH["semax_pr_CE"] = DCdf["semax_DCs_CEs"]*DCdf["maskDCBHpr_CEpos"]

    DCdf_singleBH["mP_sc_CE"] = DCdf["mP_DCs_CEs"]*DCdf["maskDCBHsc_CEpos"]
    DCdf_singleBH["mS_sc_CE"] = DCdf["mS_DCs_CEs"]*DCdf["maskDCBHsc_CEpos"]
    DCdf_singleBH["semax_sc_CE"] = DCdf["semax_DCs_CEs"]*DCdf["maskDCBHsc_CEpos"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    DCdf_tot["mP_tot"] = np.concatenate([DCdf_singleBH["mP_pr"], DCdf_singleBH["mS_sc"]])
    DCdf_tot["mS_tot"] = np.concatenate([DCdf_singleBH["mP_sc"], DCdf_singleBH["mS_pr"]])
    DCdf_tot["semax_tot"] = np.concatenate([DCdf_singleBH["semax_pr"] , DCdf_singleBH["semax_sc"] ])

    DCdf_tot["mP_tot_CE"] = np.concatenate([DCdf_singleBH["mP_pr_CE"], DCdf_singleBH["mS_sc_CE"]])
    DCdf_tot["mS_tot_CE"] = np.concatenate([DCdf_singleBH["mP_sc_CE"], DCdf_singleBH["mS_pr_CE"]])
    DCdf_tot["semax_tot_CE"] = np.concatenate([DCdf_singleBH["semax_pr_CE"] , DCdf_singleBH["semax_sc_CE"] ])
##########

    SPdf_singleBH = pd.DataFrame()
    SPdf_tot = pd.DataFrame()
    SPdf_singleBH["mP"] = SPdf["mP_innms"]*SPdf["maskSPBH"]
    SPdf_singleBH["mS"] = SPdf["mS_innms"]*SPdf["maskSPBH"]
    SPdf_singleBH["semax"] = SPdf["semax_innms"]*SPdf["maskSPBH"]

    SPdf_singleBH["mP_pr"] = SPdf["mP_innms"]*SPdf["maskSPBHpr"]
    SPdf_singleBH["mS_pr"] = SPdf["mS_innms"]*SPdf["maskSPBHpr"]
    SPdf_singleBH["semax_pr"] = SPdf["semax_innms"]*SPdf["maskSPBHpr"]

    SPdf_singleBH["mP_sc"] = SPdf["mP_innms"]*SPdf["maskSPBHsc"]
    SPdf_singleBH["mS_sc"] = SPdf["mS_innms"]*SPdf["maskSPBHsc"]
    SPdf_singleBH["semax_sc"] = SPdf["semax_innms"]*SPdf["maskSPBHsc"]

    SPdf_singleBH["mP_CE"] = SPdf["mP_ins"]*SPdf["maskSPBH_CEpos"]
    SPdf_singleBH["mS_CE"] = SPdf["mS_ins"]*SPdf["maskSPBH_CEpos"]
    SPdf_singleBH["semax_CE"] = SPdf["semax_ins"]

    SPdf_singleBH["mP_pr_CE"] = SPdf["mP_ins"]*SPdf["maskSPBHpr_CEpos"]
    SPdf_singleBH["mS_pr_CE"] = SPdf["mS_ins"]*SPdf["maskSPBHpr_CEpos"]
    SPdf_singleBH["semax_pr_CE"] = SPdf["semax_ins"]*SPdf["maskSPBHpr_CEpos"]

    SPdf_singleBH["mP_sc_CE"] = SPdf["mP_ins"]*SPdf["maskSPBHsc_CEpos"]
    SPdf_singleBH["mS_sc_CE"] = SPdf["mS_ins"]*SPdf["maskSPBHsc_CEpos"]  
    SPdf_singleBH["semax_sc_CE"] = SPdf["semax_ins"]*SPdf["maskSPBHsc_CEpos"] 

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    SPdf_tot["mP_tot"] = np.concatenate([SPdf_singleBH["mP_pr"], SPdf_singleBH["mS_sc"]])
    SPdf_tot["mS_tot"] = np.concatenate([SPdf_singleBH["mP_sc"], SPdf_singleBH["mS_pr"]])
    SPdf_tot["semax_tot"] = np.concatenate([SPdf_singleBH["semax_pr"] , SPdf_singleBH["semax_sc"] ])

    SPdf_tot["mP_tot_CE"] = np.concatenate([SPdf_singleBH["mP_pr_CE"], SPdf_singleBH["mS_sc_CE"]])
    SPdf_tot["mS_tot_CE"] = np.concatenate([SPdf_singleBH["mP_sc_CE"], SPdf_singleBH["mS_pr_CE"]])
    SPdf_tot["semax_tot_CE"] = np.concatenate([SPdf_singleBH["semax_pr_CE"] , SPdf_singleBH["semax_sc_CE"] ])
##########


    def mylog(x,y):
        return np.log10(x[y > 0]), np.log10(y[y > 0])
    
######## UPDATE 07.11.24
    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["semax_tot"],SPdf_tot["mP_tot"]), bins=20)
    df = pd.DataFrame({
    "mP_tot": np.repeat(xbins[:-1], len(ybins)-1),
    "semax_tot": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="mP_tot", columns="semax_tot", values="frequency")
    f, ax = plt.subplots(figsize=(27, 24))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["semax_tot_CE"],SPdf_tot["mP_tot_CE"]), bins=20)
    df = pd.DataFrame({
    "mP_tot_CE": np.repeat(xbins[:-1], len(ybins)-1),
    "semax_tot_CE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="mP_tot_CE", columns="semax_tot_CE", values="frequency")
    f, ax = plt.subplots(figsize=(27, 24))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Mp_CE_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["semax_tot"],SPdf_tot["mS_tot"]), bins=20)
    # df = pd.DataFrame({
    # "mS_tot": np.repeat(xbins[:-1], len(ybins)-1),
    # "semax_tot": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="mS_tot", columns="semax_tot", values="frequency")
    # f, ax = plt.subplots(figsize=(27, 24))
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

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["semax_tot_CE"],SPdf_tot["mS_tot_CE"]), bins=20)
    # df = pd.DataFrame({
    # "mS_tot_CE": np.repeat(xbins[:-1], len(ybins)-1),
    # "semax_tot_CE": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="mS_tot_CE", columns="semax_tot_CE", values="frequency")
    # f, ax = plt.subplots(figsize=(27, 24))
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
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_SP_Ms_CE_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["semax_tot"],DCdf_tot["mP_tot"]), bins=20)
    df = pd.DataFrame({
    "mP_tot": np.repeat(xbins[:-1], len(ybins)-1),
    "semax_tot": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="mP_tot", columns="semax_tot", values="frequency")
    f, ax = plt.subplots(figsize=(27, 24))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems(DCO)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["semax_tot_CE"],DCdf_tot["mP_tot_CE"]), bins=20)
    df = pd.DataFrame({
    "mP_tot_CE": np.repeat(xbins[:-1], len(ybins)-1),
    "semax_tot_CE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="mP_tot_CE", columns="semax_tot_CE", values="frequency")
    f, ax = plt.subplots(figsize=(27, 24))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCO)",  fontsize=40, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Mp_CE_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["semax_tot"],DCdf_tot["mS_tot"]), bins=20)
    # df = pd.DataFrame({
    # "mS_tot": np.repeat(xbins[:-1], len(ybins)-1),
    # "semax_tot": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="mS_tot", columns="semax_tot", values="frequency")
    # f, ax = plt.subplots(figsize=(27, 24))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCO)",  fontsize=40, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["semax_tot_CE"],DCdf_tot["mS_tot_CE"]), bins=20)
    # df = pd.DataFrame({
    # "mS_tot_CE": np.repeat(xbins[:-1], len(ybins)-1),
    # "semax_tot_CE": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="mS_tot_CE", columns="semax_tot_CE", values="frequency")
    # f, ax = plt.subplots(figsize=(27, 24))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCO)",  fontsize=40, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_" +  model + "_DC_Ms_CE_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    #####