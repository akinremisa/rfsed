"""
..module:: WaveformFittingMultiproc
        :synopsis: Calculate the Waveform Fitting for each Moho depth in parallel
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""
#%%
import numpy as np
import os 
from os.path import exists
from os import mkdir
from glob import glob
from rf import read_rf, RFStream
import matplotlib.pyplot as plt
import rffw
import math
from obspy.signal.tf_misfit import pg
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D as HL 
from sklearn.metrics import mean_squared_error
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
import multiprocessing
import itertools
#----------------------------------------------------------
def WaveformPara(rfstream, Sedthick, VpSed, VpCrust, rayp, KMoho, HMoho, gauparameter):
    """
    Calculate the Model Parameters for the Waveform Fitting

    :param rfstream: Stream of Receiver Functions
    :type rfstream: obspy Stream
    :param Sedthick: Sediment Layer Thickness
    :type Sedthick: float
    :param VpSed: P-wave velocity in the Sediment
    :type VpSed: float
    :param VpCrust: P-wave velocity in the Crust
    :type VpCrust: float
    :param rayp: Ray Parameter (s/km)
    :type rayp: float
    :param KMoho: numpy array of Vp/Vs ratios for the Moho
    :type KMoho: numpy array    
    :param HMoho: numpy array of Moho depths
    :type HMoho: numpy array
    :param gauparameter: Gaussian Parameter used in the receiver function calculation (e.g. 1.25)
    :type gauparameter: float
    :return: Dictionary of Model Parameters
    """
    VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
    SedDen=((1.6612*(VpSed)) - (0.4721 * ((VpSed)**2)) + (0.0671* ((VpSed)**3)) - (0.0043* ((VpSed)**4)) + (0.000106* ((VpSed)**5)))
    VsCrust=0.7858 - 1.2344*VpCrust + 0.7949*VpCrust**2 - 0.1238*VpCrust**3 + 0.0064*VpCrust**4
    CrustDen=((1.6612*(VpCrust)) - (0.4721 * ((VpCrust)**2)) + (0.0671* ((VpCrust)**3)) - (0.0043* ((VpCrust)**4)) + (0.000106* ((VpCrust)**5)))
    #----------------------------------------------------------
    # RF Stream 
    delta = rfstream[0].stats.delta
    b = rfstream[0].stats.starttime - rfstream[0].stats.onset
    t = (np.arange(0, len(rfstream[0].data), 1) * delta) + b
    staname = rfstream[0].stats.station
    Stacked=rfstream.select(component='R').stack()
    StackedRF=Stacked[0]
    StackedData=StackedRF.data
    Stdd=np.std(StackedData, axis=0)
    StddPlus = StackedData + Stdd
    StddMinus = StackedData - Stdd
    delay=abs(t[0])
    n=len(StackedData)
    # Waveform Fitting start from t = 0 to t = 25 seconds
    tzero=np.where(t==0)[0]
    tzero=int(tzero)  
    start=tzero
    end= 25/delta
    end=tzero + int(end)
    #----------------------------------------------------------
    parameters=[{'rfstream':rfstream, 'rayp': rayp, 'gauparameter':gauparameter, 'Sedthick':Sedthick,  'VpSed':VpSed, 'VpCrust':VpCrust, 'VsSed':VsSed, 'SedDen':SedDen, 'VsCrust':VsCrust,\
            'CrustDen':CrustDen, 'HMoho': HMoho, 'KMoho': KMoho, 'delta':delta, 'time':t, 'staname':staname, 'StackedData':StackedData, 'Stdd':Stdd, \
            'StddPlus':StddPlus, 'StddMinus':StddMinus, 'delay': delay, 'n':n, 'tstart': start, 'tend':end}]
    return parameters
#----------------------------------------------------------
def WaveformFit_multiproc(inputparams):
    """
    Calculate the Waveform Fitting for each Moho depth

    :param inputparams: List of input parameters (HMoho, OtherParams)
    :type inputparams: list
    :return: List of Waveform Fitting Results
    """
    HMoho, OtherParams = inputparams
    HMoho=HMoho[0]
    OtherParams=OtherParams[0]
    Sedthick=OtherParams['Sedthick']
    VpSed=OtherParams['VpSed']
    VsSed=OtherParams['VsSed']
    SedDen=OtherParams['SedDen']
    VpCrust=OtherParams['VpCrust']
    VsCrust=OtherParams['VsCrust']
    CrustDen=OtherParams['CrustDen']
    rayp=OtherParams['rayp']
    KMoho=OtherParams['KMoho']
    delta=OtherParams['delta']
    delay=OtherParams['delay']
    t=OtherParams['time']
    start=OtherParams['tstart']
    end=OtherParams['tend']
    n=OtherParams['n']
    StackedData=OtherParams['StackedData']
    gaussian=OtherParams['gauparameter']
    MohoSol=[]
    for i in range(len(KMoho)):
        Ktemp = KMoho[i]
        Htemp = HMoho
        #----------------------------------------------------------
        # Calculate the model parameters
        #----------------------------------------------------------
        depth = np.array([Sedthick, Htemp, 77.5])
        # print(depth)
        vp = np.array([VpSed, VpCrust, 8.045])
        vs = np.array([VsSed, (VpCrust/Ktemp), 4.49])
        rho = np.array([SedDen, CrustDen, 3.299])
        rf_synth = rffw.rffw(depth, vp, rho, vs, rayp, gaussian, delay, n, delta)
        StackedData = StackedData
        #----------------------------------------------------------
        # Goodness of Fit Evaluation
        #----------------------------------------------------------
        # general constants for goodness of fit evalaution
        dt = delta
        fmin = .5
        fmax = 10
        nf = 100
        n=n
        #Subset the Real and Synthetic Data to 0 to 25 seconds
        RealData = StackedData[start:end]
        SynthData = rf_synth[start:end]
        t=t[start:end]
        #----------------------------------------------------------
        st1a= RealData
        st2= SynthData  # st2 is the reference data
        #----------------------------------------------------------
        # Correlation Coefficient
        Correlation = np.corrcoef(RealData, SynthData)
        Correlation = Correlation[0,1]
        #----------------------------------------------------------
        # RMSE Criteria
        RMSE = math.sqrt(mean_squared_error(RealData, SynthData))
        #----------------------------------------------------------
        # Normalise Phase Goodness of Fit Criteria (PG)
        abs_PG=abs(RealData)
        abs_rf=abs(SynthData)
        PG_Eq=np.where(abs_PG == 0, abs_PG, RealData/abs_PG)
        Synth_Eq=np.where(abs_rf == 0, abs_rf, SynthData/abs_rf)
        PG=pg(PG_Eq, Synth_Eq, dt=dt, fmin=fmin, fmax=fmax, nf=nf, w0=6, st2_isref=True, a=10., k=1.)
        #----------------------------------------------------------
        MohoSol.append([Correlation, RMSE, PG, Htemp, Ktemp])
    return MohoSol
#----------------------------------------------------------
def run_waveformfitting(nproc, HMoho, ModelParams):
    """
    Run the Waveform Fitting Function in parallel for each Moho depth

    :param nproc: Number of processors to use
    :type nproc: int
    :param HMoho: numpy array of Moho depths
    :type HMoho: numpy array
    :param ModelParams: Dictionary of Model Parameters from function WaveformPara
    :type ModelParams: dict
    :return: List of Waveform Fitting Results
    """
    nproc = nproc
    HMoho = HMoho
    with multiprocessing.Pool(nproc) as pool:
        inputparams=((inputparams, ModelParams) for inputparams in itertools.product(HMoho))
        WaveformFitResults=pool.map(WaveformFit_multiproc, inputparams)
    pool.close()
    pool.join()
    return WaveformFitResults
#----------------------------------------------------------
def plotbestmodel(WaveformResults, ModelParams, wtCorr, wtRMSE, wtPG, savepath, format):
    """
    Plot the best model from the Waveform Fitting Results

    :param WaveformResults: List of Waveform Fitting Results
    :type WaveformResults: list
    :param ModelParams: Dictionary of Model Parameters
    :type ModelParams: dict
    :param wtCorr: Weighting for the Correlation Coefficient
    :type wtCorr: float
    :param wtRMSE: Weighting for the Root-mean-square Error
    :type wtRMSE: float
    :param wtPG: Weighting for the Phase Goodness of Fit
    :type wtPG: float
    :param savepath: Path to save the plot
    :type savepath: str
    :param format: Format to save the plot
    :type format: str
    """
    HSed=ModelParams[0]['Sedthick']
    VpSed=ModelParams[0]['VpSed']
    VsSed=ModelParams[0]['VsSed']
    SedDen=ModelParams[0]['SedDen']
    HMoho=ModelParams[0]['HMoho']
    KMoho=ModelParams[0]['KMoho']
    VpCrust=ModelParams[0]['VpCrust']
    VsCrust=ModelParams[0]['VsCrust']
    CrustDen=ModelParams[0]['CrustDen']
    rayp=ModelParams[0]['rayp']
    delta=ModelParams[0]['delta']
    delay=ModelParams[0]['delay']
    staname=ModelParams[0]['staname']
    gaussian=ModelParams[0]['gauparameter']
    rfstream=ModelParams[0]['rfstream']
    if wtCorr==None: wtCorr=0.4
    if wtRMSE==None: wtRMSE=0.4
    if wtPG==None: wtPG=0.2
    AllSolutions = np.zeros(((len(HMoho)*len(KMoho)), 5))
    z=0
    for i in WaveformResults:
        for j in i:
            # print(j)
            AllSolutions[z,:] = [j[0], j[1], j[2], j[3], j[4]]
            z=z+1
    #----------------------------------------------------------
    # Normalise the Correlation and RMSE and PG values 
    #----------------------------------------------------------
    Corr=AllSolutions[:,0]
    MaxCorr=max(Corr)
    RMSE=AllSolutions[:,1]
    MaxRSME=max(RMSE)
    PG=AllSolutions[:,2]
    MaxPG=max(PG)
    NormCorr=Corr/MaxCorr
    NormRMSE=RMSE/MaxRSME
    NormPG=PG/MaxPG
    GoF= wtCorr*NormCorr + wtRMSE*(1-NormRMSE) + wtPG*NormPG #Invert the RMSE values
    #-------------------------
    Htemp=AllSolutions [:,3]
    Ktemp=AllSolutions [:,4]
    #-------------------------
    AllCriteria=zip(Corr, RMSE, PG, GoF, Htemp, Ktemp)
    Solutions = np.zeros((len(GoF), 6))
    z=0
    for i in AllCriteria:
        Solutions[z,:] = [i[0], i[1], i[2], i[3], i[4], i[5]]
        z=z+1
    BestModel_Correlation = Solutions[np.argmax(Solutions[:,0]),:]
    BestModel_RMSE = Solutions[np.argmin(Solutions[:,1]),:]
    BestModel_PG = Solutions[np.argmax(Solutions[:,2]),:]
    BestModel_GoF = Solutions[np.argmax(Solutions[:,3]),:]
    print('Best Model from Correlation Coefficient: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_Correlation[4],BestModel_Correlation[5]))
    print('Best Model from RMSE: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_RMSE[4],BestModel_RMSE[5]))
    print('Best Model from PG: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_PG[4],BestModel_PG[5]))
    print('Best Model from Overall Goodness of Fit: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_GoF[4],BestModel_GoF[5]))
    labels=['Correlation Coefficient', 'Root-mean-square Error', 'Phase Goodness of Fit', 'Overall Goodness of Fit']
    model_criteria=zip(Corr, RMSE, PG, GoF)
    for i in range(len(labels)):
        suff=labels[i]
        if i==0:
            Criteria=Corr
            BestModel=BestModel_Correlation
        elif i==1:
            Criteria=RMSE
            BestModel=BestModel_RMSE
        elif i==2:
            Criteria=PG
            BestModel=BestModel_PG
        elif i==3:
            Criteria=GoF
            BestModel=BestModel_GoF
        # ----------------------------------------------------------
        # Plot the Solutions
        # ----------------------------------------------------------
        plt.figure(figsize=(18, 6)) 
        gs_corr = gridspec.GridSpec(1, 3, width_ratios=[6,8,4])
        minVpVs=min(KMoho)
        maxVpVs=max(KMoho)
        minMoho=min(HMoho)
        maxMoho=max(HMoho)
        ax1 = plt.subplot(gs_corr[0, 0])
        plt.tricontourf(Solutions[:,4], Solutions[:,5], Criteria,50, cmap ="jet")
        cb=plt.colorbar(format='%.3f', orientation="vertical")
        cb.ax.set_xlabel('S', rotation=0)
        p1,=ax1.plot(BestModel[4],BestModel[5], 'k+', color='black', mew=2, ms=10, \
        label='%s Best Model %.2f km %.2f Vp/Vs'%(staname,BestModel[4],BestModel[5]+0.001))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
            ncol=2, mode="expand", borderaxespad=0, \
            handler_map={p1:HL(numpoints=1)})
        ax1.set_ylabel('Vp/Vs')
        ax1.set_xlabel('depth (km)')
        ax1.set_ylim(minVpVs, maxVpVs)
        ax1.set_xlim(minMoho, maxMoho)
        #----------------------------------------------------------
        # Plot Stacked Data, Standard Deviation and Synthetic RF (Show time between -5 to 25 seconds)
        #----------------------------------------------------------
        delta = rfstream[0].stats.delta
        b = rfstream[0].stats.starttime - rfstream[0].stats.onset
        t = (np.arange(0, len(rfstream[0].data), 1) * delta) + b
        staname = rfstream[0].stats.station
        Stacked=rfstream.select(component='R').stack()
        StackedRF=Stacked[0]
        StackedData=StackedRF.data
        Stdd=np.std(StackedData, axis=0)
        StddPlus = StackedData + Stdd
        StddMinus = StackedData - Stdd
        n=len(StackedData)
        #----------------------------------------------------------
        H=BestModel[4]
        K=BestModel[5]
        depth = np.array([HSed, H, 77.5])
        vp = np.array([VpSed, VpCrust, 8.045])
        vs = np.array([VsSed, (VpCrust/K), 4.49])
        rho = np.array([SedDen, CrustDen, 3.299])
        rf_synth = rffw.rffw(depth, vp, rho, vs, rayp, gaussian, delay, n, delta)
        l = len(rf_synth)
        synth_time = np.arange(0, l)
        synth_time = (delta *  synth_time) - delay 
        # # Plot to visualise difference
        ax2=plt.subplot(gs_corr[0, 1])
        plt.plot(t, StackedData, color = 'black', label = 'Stacked RF')
        plt.plot(synth_time, rf_synth, '--',  color = 'black', label = 'Synthetic RF')
        plt.fill_between(t, StddPlus, StddMinus, alpha = 0.3, color = 'grey', 
                            label = 'RF Data Std' )
        ax2.set_title("%s Best Model" %suff)
        # Subset the Data to -5 to 25 seconds
        ax2.set_xlabel('time (s)')
        ax2.set_ylabel('amplitude')
        ax2.set_yticklabels([])
        ax2.set_yticks([])
        ax2.set_xlim(-5, 25)
        ax2.legend()
        #----------------------------------------------------------
        # Plot Vp and Vs with Depth
        #----------------------------------------------------------
        ax3=plt.subplot(gs_corr[0, 2])
        dep=(0, HSed, HSed, H, H, 50)
        vp=(VpSed, VpSed, VpCrust, VpCrust, 8.045, 8.045)
        vs=(VsSed, VsSed, VpCrust/K, VpCrust/K, 4.49, 4.49)
        ax3.plot(vp, dep, '--', color = 'red', label = 'Vp')
        ax3.plot(vs, dep, ':', color = 'blue', label = 'Vs')
        ax3.xaxis.set_ticks_position('bottom')
        # ax3.yaxis.tick_right()
        ax3.set_ylim(0, 50)                   
        ax3.set_xlim(1, 10)  
        ax3.invert_yaxis()
        ax3.set_ylabel('depth(km)')
        ax3.legend(fontsize=15)
        ax3.set_xlabel("v (km/s)")
        #------------------------------------------
        #    Save Plot to File
        #-------------------------------------------
        filename = "%s/%s_%s.%s"%(savepath, staname, suff, format)
        plt.savefig(filename , format=format, transparent=False, 
                    dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close("all")  
#----------------------------------------------------------