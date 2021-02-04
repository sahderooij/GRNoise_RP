import numpy as np

from kidata import io, calc, plot

import kidcalc
import scipy.integrate as integrate
from scipy import interpolate
from scipy.optimize import root
import matplotlib.pyplot as plt
import matplotlib
import warnings

def get_KIDparam(Chipnum,KIDnum,Pread,plottesc=False,
                 tescPread='max',tescrelerrthrs=.2):
    S21data = io.get_S21data(Chipnum,KIDnum,Pread)
    V = S21data[0,14]
    kbTc = S21data[0,21]*86.17

    tesc = calc.tesc(Chipnum,KIDnum,pltkaplan=plottesc,
                     Pread=tescPread,relerrthrs=tescrelerrthrs)
    return V,kbTc,tesc

#####################################################################################
##################################### Models ########################################
#####################################################################################

class Model:
    def __init__(self,V,kbTc,tesc):
        self.V = V
        self.kbTc = kbTc
        self.tesc = tesc
        #Standard constants
        self.N0 = 1.72e4
        self.t0 = 440e-3
        self.tpb = .28e-3
        self.kbTD = 37312.
        self.kb = 86.17
    
    def calc_spec(self,*args,lvlcal=1,startf=1,stopf=1e6,points=100,
                  retnumrates=False,PSDs='NN'):
        if retnumrates:
            M,B,num,rates = self.calc_MB(*args,retnumrates=True)
        else:
            M,B = self.calc_MB(*args)
        frqarr = np.logspace(np.log10(startf),np.log10(stopf),points)
        warr = 2*np.pi*frqarr*1e-6 #note the conversion to µs^-1
        sw = np.zeros(len(warr))
        for j in range(len(warr)):
            Gw = 2*np.real(np.linalg.multi_dot([
                np.linalg.inv(M+1j*warr[j]*np.eye(len(M))),
                B,
                np.linalg.inv(M.transpose()-1j*warr[j]*np.eye(len(M)))]))
#             Gw = 2/warr[j]**2*np.real(
#                 np.linalg.inv(np.eye(len(M))+M/(1j*warr[j])).dot(B))
            if PSDs == 'NN+NNt':
                sw[j] = (Gw[0,0] + Gw[1,0])
            elif PSDs == 'NN':
                sw[j] = Gw[0,0]
            elif PSDs == 'NtNt':
                sw[j] = Gw[1,1]
            elif PSDs == 'NNt':
                sw[j] = Gw[0,1]
            elif PSDs == 'NCNC':
                sw[j] = (Gw[0,0] + Gw[1,0] + Gw[0,1] + Gw[1,1])
            elif PSDs == 'N':
                sw[j] = Gw
            else:
                raise ValueError('{} is not a valid PSDs value'.format(PSDs))
        swdB = 10*np.log10(sw*1e-6*lvlcal)
        if retnumrates:
            return frqarr,swdB,num,rates
        else:
            return frqarr,swdB
    
    def calc_ltnlvl(self,Tmin,Tmax,
                    *args,lvlcal=1,points=20,
                    plotspec=False,plotnumrates=False,PSDs='NN'):
        tau = np.full(points,np.nan)
        tauerr = np.full(points,np.nan)
        lvl = np.full(points,np.nan)
        lvlerr = np.full(points,np.nan)
        Temp = np.linspace(Tmin,Tmax,points)

        cmap = matplotlib.cm.get_cmap('viridis')
        norm = matplotlib.colors.Normalize(vmin=Temp.min(),vmax=Temp.max())
        for i in range(len(Temp)):
            if plotnumrates:
                freq,swdB,nums,rates = self.calc_spec(*args,Temp[i]*self.kb,lvlcal=lvlcal,
                                                      retnumrates=True,
                                                      PSDs=PSDs)
                plt.figure('Nums')
                numcol = ['b','g','r','c','m','y','k']
                for val in nums.values():
                    plt.plot(Temp[i],val,numcol.pop(0)+'.')
                
                plt.figure('Rates')
                ratecol = ['b','g','r','c','m','y','k']
                for val in rates.values():
                    plt.plot(Temp[i],val,ratecol.pop(0)+'.')
                
            else:
                freq,swdB = self.calc_spec(*args,Temp[i]*self.kb,lvlcal=lvlcal,PSDs=PSDs)

            tau[i],tauerr[i],lvl[i],lvlerr[i] = calc.tau(freq,swdB,startf=1e0,stopf=1e5,
                                                             plot=False,retfnl=True)
            if plotspec:
                plt.figure('Spectra')
                plt.plot(freq,swdB,color=cmap(norm(Temp[i])))
        if plotspec:
            plt.figure('Spectra')
            plt.xscale('log')
            plt.ylabel('Noise level (dBc/Hz)')
            plt.xlabel('Frequency (Hz)')
            clb = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap))
            clb.ax.set_title('T (K)')
            
        if plotnumrates:
            plt.figure('Nums')
            plt.title('Numbers')
            plt.yscale('log')
            plt.xlabel('Temperature (K)')
            plt.ylabel('Number')
            plt.legend(list(nums.keys()))
            plt.figure('Rates')
            plt.title('Rates')
            plt.yscale('log')
            plt.legend(list(rates.keys()))
            plt.ylabel(r'Rate ($\mu s^{-1}$)')
            plt.xlabel('Temperature (K)')
        return Temp,tau,tauerr,lvl,lvlerr
    
    def calc_Nqpevol(self,dNqp,tStop,tInc,*args,Nqpnum='free'):
        t = np.arange(0., tStop, tInc)
        params,*rest = self.calc_params(*args)
        Nss = root(self.rateeq,
                   np.ones(self.nrRateEqs)*dNqp,
                   args=(t,params),
                   tol=1e-12,
                   jac=self.jac,method='hybr').x
        Nini = Nss.copy()
        Nini[0] += dNqp
        Nevol = integrate.odeint(self.rateeq,Nini,t,args=(params,))
        if Nqpnum == 'total':
            Nevol[:,0] += Nevol[:,1]
        return t,Nevol
    
    #plot fuctions
    def plot_ltlvl(self,T,tau,tauerr,lvl,lvlerr,ax1=None,ax2=None,color='b',fmt='',label=None):
        if ax1 is None or ax2 is None:
            fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,3))
        mask = np.logical_and(lvl/lvlerr > 2,tau/tauerr > 2)
        ax1.errorbar(T[mask]*1e3,tau[mask],yerr=tauerr[mask],color=color,fmt=fmt,label=label)
        ax1.set_yscale('log')
        ax1.set_ylabel(r'Lifetime (µs)')
        ax1.set_xlabel(r'T (mK)')
        ax2.errorbar(T[mask]*1e3,10*np.log10(lvl[mask]),
            yerr=10*np.log10((lvlerr[mask]+lvl[mask])/lvl[mask]),color=color,fmt=fmt,label=label)
        ax2.set_ylabel('FNL (dB/Hz)')
        ax2.set_xlabel(r'T (mK)')
        
#################################### Rt ##################################################
class Rt(Model):
    def __init(self,V,kbTc,tesc):
        self.nrRateEqs = 1
        super().__init__(V,kbTc,tesc)
    
    def rateeq(self,Ns,t,params):
        N = Ns
        Rstar,V,NT,NtT,Rtstar = params
        return -Rstar*(N**2-NT**2)/V-2*Rtstar*NtT*(N-NT)/V
    
    def jac(self,Ns,t,params):
        N = Ns
        Rstar,V,NT,NtT,Rtstar = params
        return -2*Rstar*N/V-2*Rtstar*NtT/V
    
    def calc_params(self,e,nrTraps,t1,xi,kbT):
        D = kidcalc.D(kbT,self.N0,self.kbTc,self.kbTD)
        c = (nrTraps/self.V/(2*self.N0*D))
        NT = kidcalc.nqp(kbT,D,self.N0)*self.V
        NtT = self.V*2*self.N0*D*c/kidcalc.f(e,kbT)
        R = (2*D/self.kbTc)**3/(4*D*self.N0*self.t0)
        Rstar = R/(1+self.tesc/self.tpb)
        Rtstar = xi*(1+e/D)/(t1*self.N0*D)
        return [Rstar,self.V,NT,NtT,Rtstar]
    
    def calc_MB(self,*args,retnumrates = False):
        params = self.calc_params(*args)
        Rstar,V,NT,NtT,Rtstar = params
        #Steady state values
        t=0 #dummy var. for rate eq.
        Nqp0 = root(
            self.rateeq,[NT],args=(t,params),
            tol=1e-12,jac=self.jac,method='hybr').x

        #M and B matrices
        M = np.array([2*Rstar*Nqp0/V+2*Rtstar*NtT/V])
        B = np.array([2*Rstar*(Nqp0**2+NT**2)/V+2*Rtstar*NtT*(Nqp0+NT)/V])

        if retnumrates:
            return M,B,{'NqpT':NT,
                        'Nqp0':Nqp0,
                        'NtT':NtT},{'R*Nqp0/2V':Rstar*Nqp0/(2*V),
                                    'Rt*NtT/2V':Rtstar*NtT/(2*V)}
        else:
            return M,B
        
####################################Trap Detrap############################################

class TrapDetrap(Model):
    def __init__(self,V,kbTc,tesc):
        self.nrRateEqs = 3
        super().__init__(V,kbTc,tesc)
    #Model equations
    def rateeq(self,Ns,t,params):
        N,Nt,Nw = Ns
        R,V,GB,Gt,Gd,Ges,NTw,eta,P,D = params
        return [-R*N**2/V + 2*GB*Nw - Gt*N + Gd*Nt + eta*P/D,
                Gt*N-Gd*Nt,
                R*N**2/(2*V) - GB*Nw - Ges*(Nw-NTw)]
    
    def jac(self,Ns,t,params):
        N,Nt,Nw = Ns
        R,V,GB,Gt,Gd,Ges,NTw,eta,P,D = params
        return [[-2*R*N/V-Gt, Gd,2*GB],
                [Gt,-Gd,0],
                [R*N/V,0,-GB-Ges]]
    
    def calc_params(self,Gt,Gd,eta,P,kbT):
        D = kidcalc.D(kbT,self.N0,self.kbTc,self.kbTD)
        NT = kidcalc.nqp(kbT,D,self.N0)*self.V
        R = (2*D/self.kbTc)**3/(4*D*self.N0*self.t0)
        Rstar = R/(1+self.tesc/self.tpb)
        Ges = 1/self.tesc
        GB = 1/self.tpb
        NTw = R*NT**2/(2*self.V*GB)
        return [R,self.V,GB,Gt,Gd,Ges,NTw,eta,P,D],NT,Rstar
    
    def calc_MB(self,Gt,Gd,eta,P,kbT,retnumrates = False):
        params,NT,Rstar = self.calc_params(Gt,Gd,eta,P,kbT)
        R,V,GB,Gt,Gd,Ges,NTw,eta,P,D = params
        #Steady state values
        t=0 #dummy var. for rate eq.
        Nqp0,Nt0,Nw0 = root(
            self.rateeq,[NT,NT/2,NTw],args=(t,params),
            tol=1e-12,jac=self.jac,method='hybr').x
        #M and B matrices
        GRstar = 2*Rstar*Nqp0/V
        M = np.array([[Gt+GRstar,-Gd],
                      [-Gt,Gd]])
        B = np.array([[2*(GRstar+Gt)*Nqp0,-1*(Gt*Nqp0+Gd*Nt0)],
                      [-1*(Gt*Nqp0+Gd*Nt0),2*Gd*Nt0]])  
        if retnumrates:
            return M,B,{'NqpT':NT,
                        'Nqp0':Nqp0,
                        'Nt0':Nt0,
                        'Nw0':Nw0},{'R*NqpT':R*NT,
                                    'R*Nqp0':R*Nqp0,
                                    'Gt':Gt,
                                    'Gd':Gd}
        else:
            return M,B