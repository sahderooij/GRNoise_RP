import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import interpolate
import warnings
import glob

import KID
from kidcalc import D, beta, cinduct, hwread, hwres, kbTeff, nqp

#################################################################
##################### GET DATA FUNCTIONS ########################    
#################################################################

def get_datafld():
    return "G:\\Thesis_data\\"

def get_S21Pread(Chipnum,KIDnum):
    datafld = get_datafld()
    S21fld = datafld + '\\'.join([Chipnum,'S21','2D'])
    return [
        int(i.split('\\')[-1].split('_')[1][:-3]) 
        for i in glob.glob(S21fld + '\\KID{}_*Tdep.csv'.format(KIDnum))]

def get_S21data(Chipnum, KIDnum, Pread):
    datafld = get_datafld()
    S21file = datafld + "\\".join(
        [
            Chipnum,
            "S21",
            "2D",
            "_".join(["KID" + str(KIDnum), str(Pread) + "dBm", "Tdep.csv"]),
        ]
    )
    S21data = np.genfromtxt(S21file, delimiter=",")[1:, :]
    return S21data

def get_Vdict(Chipnum):
    KIDlist = get_grKIDs(Chipnum)
    Volumes = {}
    for KIDnum in KIDlist:
        S21data = get_S21data(Chipnum,KIDnum,
                                       get_S21Pread(Chipnum,KIDnum)[0])
        Volumes[KIDnum] = S21data[0,14]
    return Volumes

def get_Pintdict(Chipnum):
    KIDlist = get_grKIDs(Chipnum)
    Pintdict = {}
    for KIDnum in KIDlist:
        S21Pread = np.array(get_S21Pread(Chipnum,KIDnum))
        Pintdict[KIDnum] = []
        for Pread in S21Pread:
            S21data =  get_S21data(Chipnum,KIDnum,Pread)
            Q = S21data[0,2]
            Qc = S21data[0,3]
            Qi = S21data[0,4]
            Pintdict[KIDnum].append(10*np.log10(10**(-1*Pread/10)*Q**2/Qc/np.pi))
    return Pintdict

def get_peakdata(Chipnum,KIDnum, Pread, Tbath, wvlngth, points = 3000):
    datafld = get_datafld()
    peakfile = datafld + "\\".join(
        [
            Chipnum,
            str(Tbath) + "mK",
            "_".join(
                ["KID" + str(KIDnum), str(Pread) + "dBm",
                 str(wvlngth), str(points) + "points"]
            ),
        ]
    )
    peakdata = scipy.io.loadmat(peakfile)
    peakdata_ph = peakdata["pulsemodelfo"][0]
    peakdata_amp = peakdata["pulsemodelfo_amp"][0]
    return peakdata_ph, peakdata_amp

def get_grKIDs(Chipnum):
    datafld = get_datafld() + '\\' + Chipnum + "\\NoiseTDanalyse\\"
    GRdata = scipy.io.loadmat(datafld + "TDresults")
    rdata = GRdata["TDparam"]
    return np.array([rdata['kidnr'][0][i][0,0] for i in range(rdata['kidnr'].size)])

def get_grTemp(Chipnum,KIDnum,Pread):
    datafld = get_datafld() + "\\" + Chipnum + "\\NoiseTDanalyse\\"
    GRdata = scipy.io.loadmat(datafld + "TDresults")
    rdata = GRdata["TDparam"]
    KIDlist = get_grKIDs(Chipnum)
    ind = np.where(KIDlist == KIDnum)[0][0]
    Preadar = rdata["Pread"][0, ind][:, 0]
    Tempar = rdata["Temp"][0, ind][np.where(Preadar == Pread), :]
    return Tempar[0,0][np.nonzero(Tempar[0,0])]
    
def get_grPread(Chipnum,KIDnum):
    datafld = get_datafld() + "\\" + Chipnum + "\\NoiseTDanalyse\\"
    GRdata = scipy.io.loadmat(datafld + "TDresults")
    rdata = GRdata["TDparam"]
    KIDlist = get_grKIDs(Chipnum)
    ind = np.where(KIDlist == KIDnum)[0][0]
    Preadar = rdata["Pread"][0, ind][:, 0]
    return Preadar[np.nonzero(Preadar)]
    
def get_grdata(Chipnum,KIDnum,Pread,Temp,spec = 'cross'):
    datafld = get_datafld() + "\\" + Chipnum + "\\NoiseTDanalyse\\"
    GRdata = scipy.io.loadmat(datafld + "TDresults")
    rdata = GRdata["TDparam"]
    KIDlist = get_grKIDs(Chipnum)
    ind = np.where(KIDlist == KIDnum)[0][0]
            
    Preadind = np.where(rdata["Pread"][0, ind][:, 0] == Pread)
    Tempind = np.where(rdata["Temp"][0, ind][Preadind, :][0,0] == Temp)
    freq = rdata["fmtotal"][0, ind][Preadind,Tempind][0, 0][0]
    if spec is 'cross':
        SPR = rdata["SPRrealneg"][0, ind][Preadind,Tempind][0, 0][0]
    elif spec is 'crosspos':
        SPR = rdata["SPRrealpos"][0, ind][Preadind,Tempind][0, 0][0]
    elif spec is 'crossimag':
        SPR = rdata["SPRimagneg"][0, ind][Preadind,Tempind][0, 0][0]
    elif spec is 'phase':
        SPR = rdata["SPPtotal"][0, ind][Preadind,Tempind][0, 0][0]
    elif spec is 'amp':
        SPR = rdata["SRRtotal"][0, ind][Preadind,Tempind][0, 0][0]
    else:
        raise ValueError('spec must be \'cross\', \'phase\' or \'amp\'.')
    return freq,SPR

#################################################################
##################### FILTER FUNCTIONS ##########################    
#################################################################

def del_ampNoise(freq,SPR,plot=False):
    #Delete -140 datapoints
    freq = freq[SPR!=-140]
    SPR = SPR[SPR!=-140]
    #Make it non-dB
    SPRn = 10**(SPR/10)
    #Substract amplifier noise
    if SPRn[np.logical_and(freq>3e4,freq<2e5)].size > 0 and \
        SPRn[np.logical_and(freq>3e2,freq<1e4)].size > 0:
        #sometimes amplifier noise is higher.. so check:
        if SPRn[np.logical_and(freq>3e4,freq<2e5)].mean() < \
            SPRn[np.logical_and(freq>3e2,freq<1e4)].mean():
            SPRn -= SPRn[np.logical_and(freq>3e4,freq<2e5)].max()
        else:
            SPRn -= SPRn[np.logical_and(freq>3e2,freq<1e4)].min()
    #filter positive 
    freqn = freq[SPRn>0]
    SPRn = SPRn[SPRn>0]
    #return to dB
    SPRn = 10*np.log10(SPRn)
    if plot:
        plt.figure()
        plt.plot(freq,SPR)
        plt.plot(freqn,SPRn)
        plt.xscale('log')
        plt.legend(['Input','Amp. noise filtered'])
    return freqn,SPRn
            
def del_1fNoise(freq,SPR,plot=False):
    #Delete -140 datapoints
    freq = freq[SPR!=-140]
    SPR = SPR[SPR!=-140]
    #Make it non-dB
    SPRn = 10**(SPR/10)
    SPRn -= freq**-.5*np.mean(SPRn[1:4])
    #filter positive 
    freqn = freq[SPRn>0]
    SPRn = SPRn[SPRn>0]
    #return to dB
    SPRn = 10*np.log10(SPRn)
    if plot:
        plt.figure()
        plt.plot(freq,SPR)
        plt.plot(freqn,SPRn)
        plt.xscale('log')
        plt.legend(['Input','1/f filtered'])
    return freqn,SPRn

def del_otherNoise(freq,SPR,plot=False):
    return del_1fNoise(*del_ampNoise(freq,SPR,plot=plot),plot=plot)
    
#################################################################
##################### CALCULATION FUNCTIONS #####################    
#################################################################

def ak(S21data, lbd0=0.092, N0=1.72e4, kbTD=37312.0,plot=False,reterr=False,method='df'):
    # Extract relevant data
    hw = S21data[:, 5] * 2 * np.pi * 0.6582 * 1e-9
    kbT = S21data[:, 1] * 86.17  # µeV

    # Set needed constants
    hw0 = hw[0]
    d = S21data[0, 25]
    kbTc = S21data[0,21] * 86.17
    D0 = 1.76 * kbTc

    # For D calculation:
    def integrand1(E, D):
        return 1 / np.sqrt(E ** 2 - D ** 2)
    Vsc = 1 / (integrate.quad(integrand1, D0, kbTD, args=(D0,))[0] * N0)
    
    if method is 'df':
        y = (hw - hw0) / hw0
    elif method is 'Qi':
        y = 1/S21data[:,4] - 1/S21data[0,4]
        
    # Mask the double measured temperatures, and only fit from 250 mK
    mask1 = np.zeros(len(y), dtype="bool")
    mask1[np.unique(np.round(S21data[:, 1], decimals=2),
                    return_index=True)[1]] = True
    mask = np.logical_and(mask1, (kbT >= 0.25 * 86.17))
    y = y[mask]

    x = np.zeros(len(y))
    i = 0
    s0 = cinduct(hw0, D(kbT[0], N0, Vsc, kbTD), kbT[0])
    for kbTi in kbT[mask]:
        D_0 = D(kbTi, N0, Vsc, kbTD)
        s = cinduct(hw[i], D_0, kbTi)
        if method is 'df':
            x[i] = (s[1] - s0[1]) / s0[1] * beta(lbd0, d, D_0, D0, kbTi)/4
        elif method is 'Qi':
            x[i] = (s[0] - s0[0]) / s0[1] * beta(lbd0, d, D_0, D0, kbTi)/2
        i += 1
    
    fit = curve_fit(lambda t, ak: ak * t, x, y)
    if plot:
        plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit[0]*x)
        plt.legend(['Data','Fit'])
        if method is 'df':
            plt.ylabel(r'$\delta f/f_0$')
            plt.xlabel(r'$\beta \delta \sigma_2/4\sigma_2 $')
        elif method is 'Qi':
            plt.ylabel(r'$\delta(1/Q_i)$')
            plt.xlabel(r'$\beta \delta \sigma_1/2\sigma_2 $')

    if reterr:
        return fit[0],fit[1]
    else:
        return fit[0]

def tau(freq, SPR, startf = None, stopf = None, plot=False,retfnl = False):
    #Filter non-values
    freq = freq[SPR!=-140]
    SPR = SPR[SPR!=-140]
    if startf is None:
        startf = 1e1
    if stopf is None:
        bdwth = np.logical_and(freq>3e2,freq<2e4)
        try:
            stopf = freq[bdwth][np.real(SPR[bdwth]).argmin()]
        except:
            stopf = 2e4
    
    # fitting a Lorentzian
    fitmask = np.logical_and(freq > startf, freq < stopf)
    fitfreq = freq[fitmask]
    if len(fitfreq) < 30:
        tau = np.nan
        tauerr = np.nan
        N = np.nan
        Nerr = np.nan
    else:
        try:
            fitPSD = 10**(np.real(SPR[fitmask]-SPR.max())/10) 
            #notice the normalization

            def Lorspec(f, t, N):
                SN = 4 * N * t / (1 + (2 * np.pi * f * t) ** 2)
                return SN

            fit = curve_fit(Lorspec, fitfreq, fitPSD,
                            bounds=([0, 0], [np.inf, np.inf]),
                            p0=(2e-4, 1e4))
            tau = fit[0][0]*1e6
            tauerr = np.sqrt(np.diag(fit[1]))[0]*1e6
            N = fit[0][1]*10**(np.real(SPR.max())/10)
            Nerr = np.sqrt(np.diag(fit[1]))[1]*10**(np.real(SPR.max())/10)
        except:
            tau = np.nan
            tauerr = np.nan
            N = np.nan
            Nerr = np.nan
    if plot:
        plt.figure()
        plt.plot(freq[SPR!=-140], np.real(SPR), 'o')
        if ~np.isnan(tau):
            plt.plot(fitfreq, 10*np.log10(Lorspec(fitfreq,tau*1e-6,N)))
        plt.xscale("log")
        
    if retfnl:
        return tau,tauerr,4*N*tau*1e-6,np.sqrt((4e-6*N*tauerr)**2 + (4e-6*Nerr*tau)**2)
    else:
        return tau,tauerr
    
def tau_peak(peakdata_ph,plot = False):
    t = (np.arange(len(peakdata_ph)) - 500) 
    peak = peakdata_ph
    
    fitmask = np.logical_and(t > 10, t < 1e3)
    t2 = t[fitmask]
    peak2 = peak[fitmask]
    fit = curve_fit(
        lambda x, a, b: b * np.exp(-x / a), t2, peak2, p0=(0.5e3, peak2[0])
    )[0]
    
    if plot:
        plt.figure()
        plt.plot(t, peak)
        plt.plot(t2, fit[1]*np.exp(-t2/fit[0]))
    return fit[0]

def tau_kaplan(Tmin,Tmax,tesc=.14e-3, 
               t0=.44,
               kb = 86.17,
               tpb=.28e-3,
               N0=1.72e4,
               kbTc = 1.2*86.17,
               kbTD=37312.0,):
    D0 = 1.76 * kbTc
    def integrand1(E, D):
        return 1 / np.sqrt(E ** 2 - D ** 2)
    Vsc = 1 / (integrate.quad(integrand1, D0, kbTD, args=(D0,))[0] * N0)
    T = np.linspace(Tmin,Tmax,100)
    taukaplan = np.zeros(len(T))
    for i in range(len(T)):
        D_ = D(kb*T[i]*1e-3, N0, Vsc, kbTD)
        nqp_ = nqp(kb*T[i]*1e-3, D_, N0)
        taukaplan[i] = t0*N0*kbTc**3/(4*nqp_*D_**2)*(1+tesc/tpb) 
    return T,taukaplan
    
def kbTbeff(S21data,tqpstar,
    t0=.44,
    kb=86.17,
    tpb=.28e-3,
    N0=1.72e4,
    kbTD=37312.0,
    tesc=0.14e-3,
    plot=False):
    
    kbTc = kb * S21data[0,21]
    D0 = 1.76 * kbTc  # µeV
    V = S21data[0, 14]

    def integrand1(E, D):
        return 1 / np.sqrt(E ** 2 - D ** 2)
    Vsc = 1 / (integrate.quad(integrand1, D0, kbTD, args=(D0,))[0] * N0)
    Nqp_0 = V * t0 * N0 * kbTc ** 3 / \
        (2 * D0 ** 2 * tqpstar) * 0.5 * (1 + tesc / tpb)
    
    return kbTeff(Nqp_0, N0, V, Vsc, kbTD, kbTc)
    
def _tesc(
    kbT,
    tqpstar,
    t0=.44,
    tpb=.28e-3,
    N0=1.72e4,
    kbTc=1.2 * 86.17,
    kbTD=37312.0
):
    '''Calculates the phonon escape time, based on tqp* via Kaplan. Times are in ns.'''
    D0 = 1.76 * kbTc
    def integrand1(E, D):
        return 1 / np.sqrt(E ** 2 - D ** 2)
    Vsc = 1 / (integrate.quad(integrand1, D0, kbTD, args=(D0,))[0] * N0)
    
    D_ = D(kbT, N0, Vsc, kbTD)
    nqp_ = nqp(kbT, D_, N0)
    return tpb*((4*tqpstar*nqp_*D_**2)/(t0*N0*kbTc**3)-1)

def tesc(Chipnum,KIDnum,
              minTemp=220,maxTemp=400,taunonkaplan=2e2,taures=1e1,
              pltfit=False,pltkaplan=False,
    t0=.44,
    kb=86.17,
    tpb=.28e-3,
    N0=1.72e4,
    kbTD=37312.0,
    defaulttesc=5.1e-5):
    
    Pread = get_grPread(Chipnum,KIDnum).min()
    kbTc = get_S21data(Chipnum,KIDnum,get_S21Pread(Chipnum,KIDnum)[0])[0,21]*kb
    S21Pread = np.array(get_S21Pread(Chipnum,KIDnum))
    Temp = get_grTemp(Chipnum,KIDnum,Pread)
    tescar = np.zeros(len(Temp))
    tqpstar = np.zeros(len(Temp))
    tqpstarerr = np.zeros(len(Temp))
    for i in range(len(Temp)):
        freq,SPR = get_grdata(Chipnum,KIDnum,Pread,Temp[i])
        freq,SPR = del_otherNoise(freq,SPR)
        tqpstar[i],tqpstarerr[i] = tau(freq,SPR,stopf=1e5,plot=pltfit)
        if pltfit:
            plt.title('{} KID{} -{} dBm T={} mK'.format(
                Chipnum,KIDnum,Pread,Temp[i]))
        
        if tqpstarerr[i]/tqpstar[i] > .2 or \
        (tqpstar[i] > taunonkaplan or tqpstar[i] < taures) or \
        (Temp[i] < minTemp or Temp[i] > maxTemp):
            tescar[i] = np.nan
        else:
            tescar[i] = _tesc(kb*Temp[i]*1e-3,tqpstar[i],
                             t0,tpb,N0,kbTc,kbTD)
    if tescar[~np.isnan(tescar)].size > 0:
        tesc1 = np.mean(tescar[~np.isnan(tescar)])
    else:
        tesc1 = np.nan

    if tesc1 < 0 or np.isnan(tesc1) or tesc1 > 1e-2:
        warnings.warn(
            'tesc ({}) is not valid and set to {} µs. {}, KID{}'.format(
                tesc1,defaulttesc,Chipnum,KIDnum))
        tesc1 = defaulttesc
    if pltkaplan:
        plt.figure()
        plt.errorbar(Temp,tqpstar,yerr=tqpstarerr,capsize=5.,fmt='o')
        mask = ~np.isnan(tescar)
        plt.errorbar(Temp[mask],tqpstar[mask],fmt='o')
        T,taukaplan = tau_kaplan(Temp[~np.isnan(tqpstar)].min(),
                   Temp[~np.isnan(tqpstar)].max(),tesc=tesc1,kbTc=kbTc)  
        plt.plot(T,taukaplan)
        plt.yscale('log')
        plt.xlabel('T (mK)')
        plt.ylabel(r'$\tau_{qp}^*$ (µs)')
        plt.legend(['Kaplan','Data', 'Selected Data'])
    return tesc1

def get_tescdict(Chipnum):
    tescdict = {}
    KIDlist = get_grKIDs(Chipnum)
    for KIDnum in KIDlist:
        tescdict[KIDnum] = tesc(Chipnum,KIDnum)
    return tescdict
        

#################################################################
#################### INITIALIZATION FUNCTIONS ###################    
#################################################################

def init_KID(Chipnum,KIDnum,Pread,Tbath,Teffmethod = 'GR',wvl = None,traps = False,
                t0=.44,
                kb=86.17,
                tpb=.28e-3,
                N0=1.72e4,
                kbTD=37312.0,
                lbd0=.092):
    S21Pread = np.array(get_S21Pread(Chipnum,KIDnum))
    closestPread = S21Pread[np.abs(S21Pread - Pread).argmin()]
    S21data = get_S21data(Chipnum,KIDnum,closestPread)
    if closestPread != Pread:
        warnings.warn('S21data at another Pread')
    
    Qc = S21data[0, 3]
    hw0 = S21data[0, 5]*2*np.pi*6.582e-4*1e-6
    kbT0 = kb * S21data[0, 1]
    V = S21data[0, 14]
    d = S21data[0, 25]
    kbTc = S21data[0,21] * kb
    ak1 = ak(S21data,lbd0,N0,kbTD)[0]
    
    tesc1 = tesc(Chipnum,KIDnum,t0=t0,kb=kb,tpb=tpb,N0=N0,kbTD=kbTD)
    
    if Teffmethod is 'GR':
        Temp = get_grTemp(Chipnum,KIDnum,Pread)
        taut = np.zeros(len(Temp))
        for i in range(len(Temp)):
            freq,SPR = get_grdata(Chipnum,KIDnum,Pread,Temp[i])
            taut[i] = tau(freq,SPR)[0]
        tauspl = interpolate.splrep(Temp[~np.isnan(taut)],taut[~np.isnan(taut)])
        tau1 = interpolate.splev(Tbath,tauspl)
    if Teffmethod is 'peak':
        peakdata_ph,peakdata_amp = get_peakdata(Chipnum,KIDnum,Pread,Tbath,wvl)
        tau1 = tau_peak(peakdata_ph)

    kbT = kbTbeff(S21data,tau1,t0,kb,tpb,N0,kbTD,tesc1)
        
    if traps:
        return KID.trKID(
            Qc=Qc, hw0=hw0, kbT0=kbT0, kbT=kbT, V=V, ak=ak1, d=d, kbTc = kbTc,tesc = tesc1)
    else:
        return KID.KID(
            Qc=Qc, hw0=hw0, kbT0=kbT0, kbT=kbT, V=V, ak=ak1, d=d, kbTc = kbTc,tesc = tesc1)
    
    
#################################################################
######################### PLOT FUNCTIONS ######################## 
#################################################################
def plot_spec(Chipnum,KIDlist=None,Pread='min',spec=['cross'],
              del1fNoise=False,delampNoise=False,Tmax=500,ax12=None):
    if KIDlist is None:
        KIDlist = get_grKIDs(Chipnum)
    
    if spec is 'all':
        specs = ['cross','amp','phase']
    elif type(spec) is list:
        specs = spec
    else:
        raise ValueError('Invalid Spectrum Selection')
        
    for KIDnum in KIDlist:
        allPread = get_grPread(Chipnum,KIDnum)
        if Pread is 'min':
            Preadar = np.array([allPread.max()])
        elif Pread is 'med':
            Preadarr = get_grPread(Chipnum,KIDnum)
            Preadar = np.array([Preadarr[np.abs(Preadarr.mean()-Preadarr).argmin()]])
        elif Pread is 'minmax':
            Preadar = [allPread.max(),
                       allPread.min()]
        elif Pread is 'all':
            Preadar = allPread
        else:
            raise ValueError('{} is not a valid Pread option'.format(Pread))
        if ax12 is None:
            fig,axs = plt.subplots(len(Preadar),len(specs),
                                   figsize=(6*len(specs),4*len(Preadar)),
                                   sharex=True,sharey=True)
            fig.suptitle('{}, KID{}'.format(Chipnum,KIDnum))
            if type(axs) is not np.ndarray:
                axs = np.array([[axs]])
            elif axs.ndim is not 2:
                axs = axs[np.newaxis,:]
        else:
            axs = ax12
        
        for ax1,_Pread in zip(range(len(Preadar)),Preadar):
            axs[ax1,0].set_ylabel('dBc/Hz')
            Temp = get_grTemp(Chipnum,KIDnum,_Pread)
            Temp = Temp[Temp<Tmax]
            cmap = matplotlib.cm.get_cmap('viridis')
            norm = matplotlib.colors.Normalize(Temp.min(),Temp.max())
            clb = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),
                               ax=axs[ax1,-1])
            clb.ax.set_title('T (mK)')
            for i in range(len(Temp)):
                for (ax2,spec) in zip(range(len(specs)),specs):
                    freq,SPR = get_grdata(Chipnum,KIDnum,_Pread,Temp[i],spec=spec)
                    if delampNoise:
                        freq,SPR = del_ampNoise(freq,SPR)
                    if del1fNoise:
                        freq,SPR = del_1fNoise(freq,SPR)
                    SPR[SPR==-140] = np.nan
                    axs[ax1,ax2].plot(freq,SPR,color=cmap(norm(Temp[i])))
                    axs[ax1,ax2].set_xscale('log')
                    axs[ax1,ax2].set_title(spec+ ', -{} dBm'.format(_Pread))
                    axs[-1,ax2].set_xlabel('Freq. (Hz)')
                    
                    
def plot_ltnlvl(Chipnum,KIDlist=None,pltPread='all',spec='cross',
                lvlcomp='',delampNoise=True,del1fNoise=True,relerrthrs=.3,
                pltKIDsep=True,pltthlvl=False,pltkaplan=False,pltthmfnl=False,
                ax12=None,color='Pread',fmt='-o',defaulttesc=.14e-3,showfit=False):
    def _make_fig(**kwargs):
        fig, axs = plt.subplots(1,2,figsize = (12,4))
        if color is 'Pread':
            cmap = matplotlib.cm.get_cmap('plasma')
            norm = matplotlib.colors.Normalize(-1.1*Preadar.max(),-.9*Preadar.min())
            clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title(r'$P_{read}$')
        elif color is 'Pint':
            cmap = matplotlib.cm.get_cmap('plasma')
            norm = matplotlib.colors.Normalize(np.array(Pintdict[KIDlist[k]]).min()*1.1,
                                               np.array(Pintdict[KIDlist[k]]).max()*.9)
            clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title(r'$P_{int}$')
        elif color is 'V':
            cmap = matplotlib.cm.get_cmap('cividis')
            norm = matplotlib.colors.Normalize(np.array(list(Vdict.values())).min(),
                                               np.array(list(Vdict.values())).max())
            clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title(r'$V (µm^3)$')
        elif color is 'KIDnum':
            cmap = matplotlib.cm.get_cmap('gist_rainbow')
            norm = matplotlib.colors.Normalize(np.array(KIDlist).min(),
                                               np.array(KIDlist).max())
            clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title('KID nr.')
        else:
            raise ValueError('{} is not a valid variable as color'.format(color))
        return fig,axs,cmap,norm
            
    if KIDlist is None:
        KIDlist = get_grKIDs(Chipnum)
        
    if color is 'Pint':
        Pintdict = get_Pintdict(Chipnum)
    elif color is 'V':
        Vdict = get_Vdict(Chipnum)
        
    if not pltKIDsep and (ax12 is None):
        fig,axs,cmap,norm = _make_fig()
    elif not pltKIDsep:
        axs = ax12
        
    for k in range(len(KIDlist)):
        if pltPread is 'min':
            Preadar = np.array([get_grPread(Chipnum,KIDlist[k]).max()])
        elif type(pltPread) is int:
            Preadar = np.array([np.sort(get_grPread(Chipnum,KIDlist[k]))[pltPread]])
        elif pltPread is 'med':
            Preadarr = get_grPread(Chipnum,KIDlist[k])
            Preadar = np.array([Preadarr[np.abs(Preadarr.mean()-Preadarr).argmin()]])
        elif pltPread is 'max':
            Preadar = np.array([get_grPread(Chipnum,KIDlist[k]).min()])
        elif pltPread is 'minmax':
            Preadar = np.array([get_grPread(Chipnum,KIDlist[k]).max(),
                                get_grPread(Chipnum,KIDlist[k]).min()])
        elif pltPread is 'minmedmax':
            Preadarr = get_grPread(Chipnum,KIDlist[k])
            Preadar = np.array([get_grPread(Chipnum,KIDlist[k]).max(),
                                Preadarr[np.abs(Preadarr.mean()-Preadarr).argmin()],
                                get_grPread(Chipnum,KIDlist[k]).min()])
        elif pltPread is 'all':
            Preadar = get_grPread(Chipnum,KIDlist[k])
        elif type(pltPread) is np.ndarray:
            Preadar = pltPread
        else:
            raise ValueError('{} is not a valid Pread selection'.format(pltPread))
            
        if pltKIDsep and ax12 is None:
            fig,axs,cmap,norm = _make_fig(Preadar=Preadar)
            fig.suptitle('KID{}, {}'.format(KIDlist[k],spec))
        elif pltKIDsep:
            axs=ax12

        if pltthlvl or 'tesc' in lvlcomp or pltkaplan:
            tesc_ = tesc(Chipnum,KIDlist[k],defaulttesc=defaulttesc)
    
        for Pread in Preadar:
            S21Pread = np.array(get_S21Pread(Chipnum,KIDlist[k]))
            closestPread = S21Pread[np.abs(S21Pread - Pread).argmin()]
            S21data = get_S21data(Chipnum,KIDlist[k],closestPread)
            akin = ak(S21data)
            V = S21data[0,14]
            Qspl = interpolate.splrep(S21data[:,1]*1e3,S21data[:,2],s=0)
            if spec is 'cross':
                Respspl = interpolate.splrep(
                    S21data[:,1]*1e3,np.sqrt(S21data[:,10]*S21data[:,18]),s=0)
            elif spec is 'amp':
                Respspl = interpolate.splrep(
                    S21data[:,1]*1e3,S21data[:,18],s=0)
            elif spec is 'phase':
                Respspl = interpolate.splrep(
                    S21data[:,1]*1e3,S21data[:,10],s=0)

            if lvlcomp is 'QakV':
                sqrtlvlcompspl = interpolate.splrep(
                    S21data[:,1]*1e3,
                    S21data[:,2]*akin/V,s=0)
            elif lvlcomp is 'QaksqrtV':
                sqrtlvlcompspl = interpolate.splrep(
                    S21data[:,1]*1e3,
                    S21data[:,2]*akin/np.sqrt(V),s=0)
            elif lvlcomp is 'QaksqrtVtesc':
                sqrtlvlcompspl = interpolate.splrep(
                    S21data[:,1]*1e3,
                    S21data[:,2]*akin/np.sqrt(V*\
                                              (1+tesc_/.28e-3)),s=0)
            elif lvlcomp is 'QaksqrtVtescTc':
                sqrtlvlcompspl = interpolate.splrep(
                    S21data[:,1]*1e3,
                    S21data[:,2]*akin/np.sqqrt(V*\
                                               (1+tesc_/.28e-3)*\
                                            (86.17*S21data[0,21])**3/\
                                               (S21data[0,15]/1.6e-19*1e6)**2),s=0)
            elif lvlcomp is 'Resp':            
                sqrtlvlcompspl = Respspl
            elif lvlcomp is '':
                sqrtlvlcompspl = interpolate.splrep(S21data[:,1]*1e3,np.ones(len(S21data[:,1])))
            else:
                raise ValueError('{} is an invalid compensation method'.format(
                    lvlcomp))
            
            Pint = 10*np.log10(10**(-1*Pread/10)*S21data[0,2]**2/S21data[0,3]/np.pi)
            Temp = np.trim_zeros(get_grTemp(Chipnum,KIDlist[k],Pread))
            taut = np.zeros((len(Temp)))
            tauterr = np.zeros((len(Temp)))
            lvl = np.zeros((len(Temp)))
            lvlerr = np.zeros((len(Temp)))
            for i in range(len(Temp)):
                freq,SPR = get_grdata(Chipnum,KIDlist[k],Pread,Temp[i],spec)
                if delampNoise:
                    freq,SPR = del_ampNoise(freq,SPR)
                if del1fNoise:
                    freq,SPR = del_1fNoise(freq,SPR)
                taut[i],tauterr[i],lvl[i],lvlerr[i] = \
                    tau(freq,SPR,plot=showfit,retfnl = True)
                if showfit:
                    plt.title('{}, KID{}, -{} dBm, T={}, {},\n relerr={}'.format(
                        Chipnum,KIDlist[k],Pread,Temp[i],spec,tauterr[i]/taut[i]))
                    
                Q = interpolate.splev(Temp[i],Qspl)
                                
                lvl[i] = lvl[i]/interpolate.splev(Temp[i],sqrtlvlcompspl)**2
                lvlerr[i] = lvlerr[i]/interpolate.splev(Temp[i],sqrtlvlcompspl)**2

            #Deleting bad fits and plotting:
            mask = ~np.isnan(taut)
            mask[mask] = tauterr[mask]/taut[mask] <= relerrthrs
            
            if color is 'Pread':
                clr = cmap(norm(-1*Pread))
            elif color is 'Pint':
                clr = cmap(norm(Pint))
            elif color is 'V':
                clr = cmap(norm(Vdict[KIDlist[k]]))
            elif color is 'KIDnum':
                clr = cmap(norm(KIDlist[k]))
            else:
                clr = color
            
            axs[0].errorbar(Temp[mask],taut[mask],
                                     yerr = tauterr[mask],fmt = fmt,capsize = 5.,
                              color=clr)
            axs[1].errorbar(Temp[mask],10*np.log10(lvl[mask]),
                                     yerr = 10*np.log10((lvlerr[mask]+lvl[mask])/lvl[mask]),
                                     fmt = fmt, capsize = 5.,color=clr)
            if pltthlvl:
                Ttemp = np.linspace(50,400,100)
                explvl = interpolate.splev(Ttemp,Respspl)**2
                explvl *= .44e-6*V*1.72e4*(86.17*S21data[0,21])**3/\
                ((S21data[0,15]/1.602e-19*1e6)**2)*(1+tesc_/.28e-3)
                explvl /= interpolate.splev(Ttemp,sqrtlvlcompspl)**2
                thlvlplot, = axs[1].plot(Ttemp,10*np.log10(explvl),color=clr,linestyle='--')
                axs[1].legend((thlvlplot,),(r'FNL from responsivity',))
            if pltkaplan:
                T,taukaplan = tau_kaplan(Temp[mask].min(),Temp[mask].max(),
                                         tesc=tesc_,kbTc=86.17*S21data[0,21])
                kaplanfit, = axs[0].plot(T,taukaplan,color=clr,linestyle='--',linewidth=3.)
                axs[0].legend((kaplanfit,),('Kaplan',))
            if pltthmfnl:
                try:
                    tauspl = interpolate.splrep(Temp[mask],taut[mask],s=0)
                    T = np.linspace(Temp[mask].min(),Temp[mask].max(),100)
                    Nqp = np.zeros(len(T))
                    for i in range(len(T)):
                        Nqp[i] = V*nqp(T[i]*86.17*1e-3,S21data[0,15]/1.602e-19*1e6,1.72e4)
                    thmfnl = 4*interpolate.splev(T,tauspl)*1e-6*\
                        Nqp*interpolate.splev(T,Respspl)**2
                    thmfnl /= interpolate.splev(T,sqrtlvlcompspl)**2
                    thmfnlplot, = axs[1].plot(T,10*np.log10(thmfnl),color=clr,
                                              linestyle='--',linewidth=3.)
                    axs[1].legend((thmfnlplot,),('Thermal FNL \n with measured $\\tau_{qp}^*$',))
                except:
                    warnings.warn('Could not make Thermal FNL, {},KID{},-{} dBm,{}'.format(
                    Chipnum,KIDlist[k],Pread,spec))
                    
        axs[0].set_title('Lifetime')
#         axs[0].set_ylim(8e0,4e3)
        axs[0].set_yscale('log')
        axs[1].set_title('Flat Noise Level')
#         axs[1].set_ylim(-160,-110)
        for i in range(2):
            axs[i].set_xlabel('T (mK)')
#             axs[i].set_xlim(Temp.min()-10,Temp.max()+10)
        axs[0].set_ylabel(r'$\tau$ (µs)')
        if lvlcomp is 'QakV':
            axs[1].set_ylabel(r'dBc/Hz-20$log_{10}$(Q$\alpha_k$/V)')
        elif lvlcomp is 'Resp':
            axs[1].set_ylabel(r'dBc/Hz-20$log_{10}$($Resp.^2)$')
        elif lvlcomp is 'QaksqrtV':
            axs[1].set_ylabel(r'dBc/Hz-20$log_{10}(Q\alpha_k/\sqrt{V})$')
        elif lvlcomp is 'QaksqrtVtesc':
            axs[1].set_ylabel(r'dBc/Hz-20$log_{10}(Q\alpha_k\sqrt{(1+\tau_{esc}/\tau_{pb})/V})$')
        elif lvlcomp is 'QaksqrtVtescTc':
            axs[1].set_ylabel(
                r'dBc/Hz-20$log_{10}(Q\alpha_k\sqrt{(1+\tau_{esc}/\tau_{pb})(k_BT_c)^3/(V\Delta^2)})$')
        else:
            axs[1].set_ylabel(r'dBc/Hz')
    
def plot_rejspec(Chipnum,KIDnum,sigma,Trange = (0,400),sepT = False, spec='SPR'):
    dfld = get_datafld()
    for j in range(len(sigma)):
        #Load matfile and extract Pread en temperature
        matfl = scipy.io.loadmat(
            dfld + '\\'+ Chipnum + '\\Noise_vs_T' + '\\TDresults_{}'.format(sigma[j])
        )['TDparam']
        if Chipnum is 'LT139':
            if KIDnum == 1:
                ind = 0
            elif KIDnum == 6:
                ind = 2
            elif KIDnum == 7:
                ind = 3
            else:
                raise ValueError('KID{} is not in LT139'.format(KIDnum))
        else:
            ind = KIDnum - 1 
        Preadar = matfl['Pread'][0,ind][:,1]
        Pread = Preadar[-1]
        Pind = np.where(Preadar == Pread)
        Tar = matfl['Temp'][0,ind][Pind,:][0,0]
        Tmin,Tmax = Trange
        Tar=Tar[np.logical_and(Tar>Tmin,Tar<Tmax)]
        Tar=np.flip(Tar) #To plot the lowest temperatures last        

        #Plot spectrum from each temperature
        if not sepT:
            plt.figure(figsize=(10,4))
            plt.subplot(1,2,1)
            cmap = matplotlib.cm.get_cmap('viridis')
            norm = matplotlib.colors.Normalize(Tar.min(),Tar.max())
            percrej = np.zeros(len(Tar))
        for i in range(len(Tar)):
            Tind = np.where(np.flip(Tar) == Tar[i])
            freq = matfl['fmtotal'][0,ind][Pind,Tind][0,0][0]
            if spec is 'SPR':
                SPR = matfl['SPRrealneg'][0,ind][Pind,Tind][0,0][0]
            elif spec is 'SPP':
                SPR = matfl['SPPtotal'][0,ind][Pind,Tind][0,0][0]
            elif spec is 'SRR':
                SPR = matfl['SRRtotal'][0,ind][Pind,Tind][0,0][0]
            else:
                raise ValueError('Choose SPR, SPP or SRR as spec')
            SPR[SPR==-140] = np.nan #delete -104 dBm, because that is minimum.
            if sepT:
                plt.figure(np.floor(i/2),figsize=(10,4))
                plt.subplot(1,2,i%2+1)
                plt.title('KID{}, -{} dBm, {} mK'.format(KIDnum,Pread,Tar[i]))
                plt.plot(freq,SPR)
                plt.legend(sigma,title=r'Thrshld ($\sigma$)')
            else:
                plt.plot(freq,SPR,color=cmap(norm(Tar[i])))
                percrej[i] = 100*matfl['nrrejectedmed'][0,ind][Pind,Tind][0,0]/32
            plt.xscale('log')
        if not sepT:
            plt.title(r'KID{}, $P_{{read}}=$-{} dBm, {}$\sigma$'.format(
                KIDnum,Pread,sigma[j]))
            plt.ylim(-130,-50)
            clb = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title('T (mK)')
            plt.subplot(1,2,2)
            plt.title('Rejection Percentage')
            plt.plot(Tar,percrej)
            plt.xlabel('T (mK)')
            plt.ylabel('%')
            plt.ylim(0,100)
            plt.tight_layout()
            

