"""
.. module:: WaveformFitting
        :synopsis: Waveform Fitting for Moho depth and Vp/Vs ratio
.. moduleauthor:: Stephen Akinremi <s.akinremi@utwente.nl> Islam Fadel <i.e.a.m.fadel@utwente.nl> (October 2023)
"""

import numpy as np
from rfsed.SynRF.FwdRF import SynRF
import math
from obspy.signal.tf_misfit import pg
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D as HL 
from sklearn.metrics import mean_squared_error
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms

def WaveformFitting(rfstream, HSed, VpSed, VpCrust, rayp, KMoho, HMoho, gaussian, wtCorr, wtRMSE, wtPG, savepath, format):
    """
    This function performs waveform fitting for the Moho depth and Vp/Vs ratio
    Returns the best model based on the Correlation Coefficient(CC), Root-mean-squared error(RMSE), 
    Phase-goodness-of-fit(PG) and Overall Goodness-of-fit that combines CC, RMSE, and PG
    Weights of CC, RMSE and PG (wtCorr,wtRMSE, wtPG) are given as input and can be modified before plotting
    At all times, wtCorr+wtRMSE+wtPG=1.0

    :param rfstream: Obspy RFStream object containing the receiver functions
    :type rfstream: Obspy RFStream object
    :param HSed: Sediment thickness in km
    :type HSed: float
    :param VpSed: Sediment Vp in km/s
    :param VpSed: float
    :param VpCrust: Crustal Vp in km/s
    :param VpCrust: float
    :param rayp: Ray parameter in s/km
    :param rayp: float
    :param KMoho: Vp/Vs ratio range for Moho
    :type KMoho: numpy array
    :param HMoho: Moho depth range in km
    :type HMoho: numpy array
    :param gaussian: Gaussian filter parameter (a)
    :type gaussian: float
    :param wtCorr: Weight of Correlation Coefficient
    :type wtCorr: float
    :param wtRMSE: Weight of Root-mean-squared error
    :type wtRMSE: float
    :param wtPG: Weight of Phase Goodness of Fit
    :type wtPG: float
    :param savepath: Path to save the plots
    :type savepath: str
    :param format: Format of the plots
    :type format: str
    :return: Dictionary containing the Waveform Fitting Results
    """
    VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
    SedDen=((1.6612*(VpSed)) - (0.4721 * ((VpSed)**2)) + (0.0671* ((VpSed)**3)) - (0.0043* ((VpSed)**4)) + (0.000106* ((VpSed)**5)))
    VsCrust=0.7858 - 1.2344*VpCrust + 0.7949*VpCrust**2 - 0.1238*VpCrust**3 + 0.0064*VpCrust**4
    CrustDen=((1.6612*(VpCrust)) - (0.4721 * ((VpCrust)**2)) + (0.0671* ((VpCrust)**3)) - (0.0043* ((VpCrust)**4)) + (0.000106* ((VpCrust)**5)))
    #----------------------------------------------------------
    # RF Stream 
    delta = rfstream[0].stats.delta
    # print('Data delta is', delta)
    b = rfstream[0].stats.starttime - rfstream[0].stats.onset
    t = (np.arange(0, len(rfstream[0].data), 1) * delta) + b
    staname = rfstream[0].stats.station
    Stacked=rfstream.select(component='R').stack()
    StackedRF=Stacked[0]
    StackedData=StackedRF.data
    Stdd=np.std(StackedData, axis=0)
    StddPlus = StackedData + Stdd
    StddMinus = StackedData - Stdd
    delay = abs(t[0])
    n=len(StackedData)
    # Waveform Fitting start from t = 0 to t = 25 seconds
    tzero=np.where(t==0)[0]
    tzero=int(tzero)  
    start=tzero
    end= 25/delta
    end=tzero + int(end)
    #----------------------------------------------------------
    MohoSol=[]
    z=0
    #----------------------------------------------------------
    for i in range(len(KMoho)):
        Ktemp = KMoho[i]
        for j in range(len(HMoho)): 
            Htemp = HMoho[j]
            #----------------------------------------------------------
            # Calculate the model parameters
            #----------------------------------------------------------
            depth = np.array([HSed, Htemp, 77.5])
            vp = np.array([VpSed, VpCrust, 8.045])
            vs = np.array([VsSed, (VpCrust/Ktemp), 4.49])
            rho = np.array([SedDen, CrustDen, 3.299])
            Synth=SynRF(depth, vp, vs, rho, rayp, dt=delta, npts=n, ipha=1)
            Synth.run_fwd()
            Synth.filter(freqmin=0.05, freqmax=1.25, order=2, zerophase=True)
            rf_synth=Synth.run_deconvolution(pre_filt=[0.05, 1.25], shift=delay, gaussian=gaussian)
            rf_synth=(rf_synth[0]).data
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
    AllSolutions = np.zeros((len(MohoSol), 5))
    for i in MohoSol:
        AllSolutions[z,:] = [i[0], i[1], i[2], i[3], i[4]]
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
    # #-------------------------
    AllCriteria=zip(Corr, RMSE, PG, GoF, Htemp, Ktemp)
    WaveformSolutions = np.zeros((len(GoF), 6))
    z=0
    for i in AllCriteria:
        WaveformSolutions[z,:] = [i[0], i[1], i[2], i[3], i[4], i[5]]
        z=z+1
    BestModel_Correlation = WaveformSolutions[np.argmax(WaveformSolutions[:,0]),:]
    BestModel_RMSE = WaveformSolutions[np.argmin(WaveformSolutions[:,1]),:]
    BestModel_PG = WaveformSolutions[np.argmax(WaveformSolutions[:,2]),:]
    BestModel_GoF = WaveformSolutions[np.argmax(WaveformSolutions[:,3]),:]
    print('Best Model from Correlation Coefficient: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_Correlation[4],BestModel_Correlation[5]))
    print('Best Model from RMSE: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_RMSE[4],BestModel_RMSE[5]))
    print('Best Model from PG: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_PG[4],BestModel_PG[5]))
    print('Best Model from Overall Goodness of Fit: Moho at %.2f km with Vp/Vs of %.2f'%(BestModel_GoF[4],BestModel_GoF[5]))
    WaveformFitting={'WaveformSolutions':WaveformSolutions, 'rfstream':rfstream, 'HSed':HSed, 'VpSed':VpSed, 'VpCrust':VpCrust, 
                     'rayp':rayp, 'KMoho':KMoho, 'HMoho':HMoho, 'gaussian':gaussian, 'wtCorr':wtCorr, 'wtRMSE':wtRMSE, 'wtPG':wtPG, 
                     'savepath':savepath, 'format':format}
    return WaveformFitting

def PlotWaveformFitting(WaveformFitting, wtCorr, wtRMSE, wtPG, savepath, format):
    """
    This function plots the waveform fitting results for the Moho depth and Vp/Vs ratio
    Plots the Waveform Fitting Results for CC, RMSE, PG and GoF
    Weights of CC, RMSE and PG are given as input and can be modified before plotting
    At all times, wtCorr+wtRMSE+wtPG=1.0

    :param WaveformFitting: Dictionary containing the Waveform Fitting Results from the function WaveformFitting
    :type WaveformFitting: dict
    :param wtCorr: Weight of Correlation Coefficient
    :type wtCorr: float
    :param wtRMSE: Weight of Root-mean-squared error
    :type wtRMSE: float
    :param wtPG: Weight of Phase Goodness of Fit
    :type wtPG: float
    :param savepath: Path to save the plots
    :type savepath: str
    :param format: Format of the plots
    :type format: str
    """
    WaveformSolutions=WaveformFitting['WaveformSolutions']
    rfstream=WaveformFitting['rfstream']
    staname = rfstream[0].stats.station
    HSed=WaveformFitting['HSed']
    VpSed=WaveformFitting['VpSed']
    VpCrust=WaveformFitting['VpCrust']
    rayp=WaveformFitting['rayp']
    KMoho=WaveformFitting['KMoho']
    HMoho=WaveformFitting['HMoho']
    gaussian=WaveformFitting['gaussian']
    BestModel_Correlation = WaveformSolutions[np.argmax(WaveformSolutions[:,0]),:]
    BestModel_RMSE = WaveformSolutions[np.argmin(WaveformSolutions[:,1]),:]
    BestModel_PG = WaveformSolutions[np.argmax(WaveformSolutions[:,2]),:]
    BestModel_GoF = WaveformSolutions[np.argmax(WaveformSolutions[:,3]),:]
    Corr = WaveformSolutions[:,0]
    RMSE = WaveformSolutions[:,1]
    PG = WaveformSolutions[:,2]
    GoF = WaveformSolutions[:,3]
    #----------------------------------------------------------
    labels=['Correlation Coefficient', 'Root-mean-square Error', 'Phase Goodness of Fit', 'Overall Goodness of Fit']
    model_criteria=zip(Corr, RMSE, PG, GoF)
    print(model_criteria)
    for i in model_criteria:
        print(i)
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
        minMoho=HMoho[0]
        maxMoho=HMoho[-1]
        minVpVs=KMoho[0]
        maxVpVs=KMoho[-1]
        # ----------------------------------------------------------
        # Plot the Solutions
        # ----------------------------------------------------------
        plt.figure(figsize=(18, 6)) 
        gs_corr = gridspec.GridSpec(1, 3, width_ratios=[6,8,4])
        ax1 = plt.subplot(gs_corr[0, 0])
        plt.tricontourf(WaveformSolutions[:,4], WaveformSolutions[:,5], Criteria,50, cmap ="jet")
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
        # Plot Stacked Data, Standard Deviation and Synthetic RF )(Show time between -5 to 25 seconds)
        #----------------------------------------------------------
        delta = rfstream[0].stats.delta
        b = rfstream[0].stats.starttime - rfstream[0].stats.onset
        t = (np.arange(0, len(rfstream[0].data), 1) * delta) + b
        delay=abs(t[0]) 
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
        # Calculate sediment and crustal density using Nafe-Drake curve (Brocher Theorem) valid for vp between 1.5 and 8.5 km/sec
        VsSed=0.7858 - 1.2344*VpSed + 0.7949*VpSed**2 - 0.1238*VpSed**3 + 0.0064*VpSed**4
        SedDen = ((1.6612*(VpSed)) - (0.4721 * ((VpSed)**2)) + (0.0671* ((VpSed)**3)) - (0.0043* ((VpSed)**4)) + (0.000106* ((VpSed)**5)))
        CrustDen = ((1.6612*(VpCrust)) - (0.4721 * ((VpCrust)**2)) + (0.0671* ((VpCrust)**3)) - (0.0043* ((VpCrust)**4)) + (0.000106* ((VpCrust)**5)))
        depth = np.array([HSed, H, 77.5])
        vp = np.array([VpSed, VpCrust, 8.045])
        vs = np.array([VsSed, (VpCrust/K), 4.49])
        rho = np.array([SedDen, CrustDen, 3.299])
        Synth=SynRF(depth, vp, vs, rho, rayp, dt=delta, npts=n, ipha=1)
        Synth.run_fwd()
        Synth.filter(freqmin=0.05, freqmax=1.25, order=2, zerophase=True)
        rf_synth=Synth.run_deconvolution(pre_filt=[0.05, 1.25], shift=delay, gaussian=gaussian)
        rf_synth=(rf_synth[0]).data
        l = len(rf_synth)
        synth_time = np.arange(0, l)
        synth_time = (delta *  synth_time) - delay 
        # # Plot to visualise difference
        ax2=plt.subplot(gs_corr[0, 1])
        plt.plot(t, StackedData, color = 'black', label = 'Stacked RF')
        plt.plot(synth_time, rf_synth, '--',  color = 'black', label = 'Synthetic RF')
        plt.fill_between(t, StddPlus, StddMinus, alpha = 0.3, color = 'grey', 
                         label = 'RF Data Std' )
        # p1,=ax2.plot(BestModel[4],BestModel[5], 'k+', color='black', mew=2, ms=10, \
        # label='%s Best Model'%(suff))
        ax2.set_title("%s Best Model" %suff)
        # Subset the Data to -5 to 25 seconds
        ax2.set_xlabel('time (s)')
        ax2.set_ylabel('amplitude')
        ax2.set_yticklabels([])
        ax2.set_yticks([])
        ax2.set_xlim(-5, 25)
        ax2.legend()
        # #----------------------------------------------------------
        # # Plot Vp and Vs with Depth
        # #----------------------------------------------------------
        ax3=plt.subplot(gs_corr[0, 2])
        dep=(0, HSed, HSed, H, H, 50)
        vp=(VpSed, VpSed, VpCrust, VpCrust, 8.045, 8.045)
        vs=(VsSed, VsSed, VpCrust/K, VpCrust/K, 4.49, 4.49)
        ax3.plot(vp, dep, '--', color = 'red', label = 'Vp')
        ax3.plot(vs, dep, ':', color = 'blue', label = 'Vs')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.set_ylim(0, 50)                   
        ax3.set_xlim(1, 10)  
        ax3.invert_yaxis()
        ax3.set_ylabel('depth(km)')
        ax3.legend(fontsize=15)
        ax3.set_xlabel("v (km/s)")#, fontsize=20)
        #------------------------------------------
        #    Save Plot to File
        #-------------------------------------------
        filename = "%s/%s_%s.%s"%(savepath, staname, suff, format)
        plt.savefig(filename , format=format, transparent=False, 
                    dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.show()
        plt.close("all")  
#----------------------------------------------------------