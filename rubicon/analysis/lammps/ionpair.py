# coding: utf-8 

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import numpy as np
from scipy.integrate import cumtrapz
from rubicon.analysis.lammps._md_analyzer import ipcorr
from rubicon.analysis.lammps._md_analyzer import calcdistances
from scipy.optimize import curve_fit
import copy
import warnings
from six.moves import range

__author__ = "mhumbert"


class ionpair:
    """
        Calculates the ionpair lifetimes for a lammps simulation.
        
        Takes in the center of mass coordinates for all molecules and timesteps 
        from calcCOM as well as other properties already parsed. 
        
        The correlation function is fit to a single, double, triple, quadrouple and
        quintuple exponential, to find the best fit without overfitting the data. 
        
        output includes the correlation function, the ionpair lifetimes from the 
        different fits as well as the r-squared values for the different fits
        
    """
    
    def runionpair(self, comx, comy, comz, Lx, Ly, Lz, moltypel, moltype, 
                   tsjump, dt, output, skipframes):
        output['Ion_Pair_Lifetime'] = {}
        output['Ion_Pair_Lifetime']['Units'] = 'picoseconds'
        output['Ion_Pair_Lifetime'][
            'Explanation'] = 'The Ion Pair Lifetime correlation function is fit to a single exponential, a double exponential up to 5 exponentials. The results shown are the result of these successive fittings'
        (closest, begin, end, C)=self.init(len(comx[0]), moltypel, len(comx),
                                           moltype)
        for step in range(0,len(comx)):
            r = self.calcdistance(comx[step], comy[step], comz[step], Lx, Ly, Lz)
            self.findclosest(r, closest, begin, end, step)
        
        
        correlation = self.correlation(closest, moltype, moltypel, skipframes)
        time = []
        for i in range(0,len(correlation)):
            time.append(float(i*tsjump*dt/1000))
        begin = int(1000/dt/tsjump)
        for i in range(0,len(moltypel)):
            for j in range(0,len(moltypel)):
                if i != j:
                    y = []
                    end = 0
                    for k in range(0,len(correlation)):
                        y.append(float(correlation[k][i][j]))
                        if correlation[k][i][j] <= 0.04 and end ==0:
                            end = k
                    if end==0:
                        end = len(y)
                    (IPL,r2) = self.curvefit(y,time,begin,len(y))
                    output['Ion_Pair_Lifetime']['{0} around {1}'.format(
                                                  moltypel[j],moltypel[i])]=IPL
                    output['Ion_Pair_Lifetime']['{0} around {1} r2'.format(
                                                   moltypel[j],moltypel[i])]=r2
                    output['Ion_Pair_Lifetime']['{0} around {1} correlation'.format(
                                         moltypel[j],moltypel[i])]=copy.deepcopy(y)                   
        output['Ion_Pair_Lifetime']['Correlation_Time']=copy.deepcopy(time)
    def init(self, nummol, moltypel, numtimesteps, moltype):
        closest = np.zeros((numtimesteps, nummol, len(moltypel)))
        C = np.zeros((len(moltypel), len(moltypel)))
        begin = [0]
        end = []
        for i in range(1, len(moltype)):
            if moltype[i]!= moltype[i-1]:
                end.append(i)
                begin.append(i)
                
        end.append(len(moltype))
        return (closest, begin, end, C)
        
    def calcdistance(self, comx, comy, comz, Lx, Ly, Lz):
        r = calcdistances(len(comx), comx, comy, comz, Lx, Ly, Lz)
        return r
        
    def findclosest(self, r, closest, begin, end, timestep):
        for i in range(0, len(r)):
            for j in range(0, len(begin)):
                distance = 10000
                for k in range(begin[j], end[j]):
                    if r[i][k]<distance:
                        distance = r[i][k]
                        closest[timestep][i][j] = k
    
    def correlation(self, closest, moltype, moltypel, skipframes):
        correlation = ipcorr(closest, skipframes, len(closest),
                                    len(closest[0]), len(closest[0][0]),
                                    (len(closest)-skipframes)/2, moltype)
        return correlation
        
    def cumulativeintegral(self, correlation, time):
        cumintegral = np.zeros((len(correlation), len(correlation[0]),
                                len(correlation[0][0])))
        for i in range(0, len(correlation[0])):
            for j in range(0, len(correlation[0][0])):
                y = []
                for k in range(0, len(correlation)):
                    y.append(correlation[k][i][j])
                cumint = cumtrapz(y, time)
                for k in range(0, len(cumint)):
                    cumintegral[k][i][j] = cumint[k]
        return cumintegral
    
    def exponential1(self, x, A1, B1):
        return A1*np.exp(-x/B1)
        
    def exponential2(self, x, A1, A2, B1, B2):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2)
        
    def exponential3(self, x, A1, A2, A3, B1, B2, B3):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3)
        
    def exponential4(self, x, A1, A2, A3, A4, B1, B2, B3, B4):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3) + (
               A4*np.exp(-x/B4))
        
    def exponential5(self, x, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3) + (
               A4*np.exp(-x/B4) + A5*np.exp(-x/B5))
        
    def curvefit(self, correlation, time, begin, end):
        funlist = [self.exponential1, self.exponential2, self.exponential3,
                   self.exponential4, self.exponential5]
        IPL = []
        r2 = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for fun in funlist:
                try:
                    popt, pcov = curve_fit(fun, np.array(time[begin:end]), 
                                           np.array(correlation[begin:end]), maxfev=1000000)
                except:
                    print("Fitting for Ion Pair Lifetime did not converge")
                    r2.append('fitting did not converge')
                    IPL.append('fitting did not converge')
                    continue
                fit = []
                for i in time:
                    fit.append(fun(i, *popt))
                yave = np.average(correlation[begin:end])
                SStot = 0
                SSres = 0
                for l in range(begin,end):
                    SStot += (correlation[l]-yave)**2
                    SSres += (correlation[l]-fit[l])**2
                r2.append(1-SSres/SStot)
                IPL.append(0)
                for i in range(0,int(len(popt)/2)):
                    IPL[-1] += popt[i]*popt[i+len(popt)/2]
            
            #plt.plot(time,correlation)
            #plt.plot(time,fit)
            #plt.show()
        return (IPL, r2)
