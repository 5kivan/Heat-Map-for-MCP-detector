# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:38:39 2020

@author: Ivan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import peakutils
from scipy.optimize import curve_fit
from scipy import special
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

#read the data and take only the meaningfull part
file_name=input('Please enter the file name of the dat format >')
data = pd.read_csv(file_name + '.dat',sep='	', decimal=',', float_precision='high', header=None)
data.columns=["energy", "counts"]
data.energy=data.energy
data[data<0]=0

# parameters of the peak estimator and loop of estimation

while True:
    threshold = float(input('The threshold for the peak estimation \n number between 0 and 1 > '))
    minimum_distanse=float(input('Minimum distance between peaks  > '))

#find the peaks
    cb = np.array(data.counts)
    peaks = peakutils.indexes(cb, thres=threshold, min_dist=minimum_distanse)
    peaks= pd.DataFrame(data, index=peaks)

#plot the data
    fig, ax = plt.subplots()
    print('Energy range to plot the graph')
    x_min=float(input('x_min > '))
    x_max=float(input('x_max > '))
    print('Counts range to plot the graph')
    y_min=float(input('y_min > '))
    y_max=float(input('y_max > '))

    ax.tick_params(axis='both',direction='out', top=True,right=True)
    ax.set_title('$^{48}Ca+^{242}Pu$', fontsize=20)
    plt.step(data.energy, data.counts,c='black',alpha=0.7)
    plt.scatter(peaks.energy,peaks.counts,c='red',marker='*',alpha=1)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_xlabel('Energy, keV', fontsize=18)
    plt.ylim(y_min, y_max)
    plt.xlim(x_min,x_max )
    plt.legend(('Data', 'Peaks'), loc='upper left')
    while True:
        isotopes = input("add the name of the isope on plot? y/n > ")
        if isotopes =="y":
           name=input("name of the isotope > ")
           x_pos= float(input("x position on the plot > "))
           y_pos=float(input("y position on the plot > "))
           ener=input("alpha energy > ")
           plt.text(x_pos, y_pos, name + '\n $E_{\u03B1}=$' + ener )
        if isotopes == "n":
           break
    saving=input('Save the plot? y/n > ')
    if saving =="y":
        plt.tight_layout()
        plt.savefig(file_name + ".pdf")
        plt.show()
    else:
        plt.show()
    cont = input("Repeat the peak finding and ploting? y/n > ")
    if cont == "n":
           break

print('Identified peaks \n')
peaks=peaks.reset_index()
peaks=peaks.drop(['index'], axis=1)
print(peaks)

#peak fitting

def line_gaus(x, amplitude, mean, stddev,slope,incr):
    return amplitude*1/(stddev*np.sqrt(2*np.pi)) * np.exp(-0.5*((x - mean)/stddev)**2)+slope*x+incr

def gaus(x, amplitude, mean, stddev):
        return amplitude*1/(stddev*np.sqrt(2*np.pi)) * np.exp(-0.5*((x - mean)/stddev)**2)

def alpha_gaus(x, amplitude, mean, stddev,tau):
        return amplitude/(2*tau) * np.exp((x - mean)/tau+0.5*(stddev/tau)**2)*special.erfc((1/np.sqrt(2))*((x - mean)/stddev+(stddev/tau)))

while True:
    x=int(input('The serial number of the peak to be fitted > '))
    a=float(input('How much of the data to fit? > '))
    left_gauss_bound = peaks.energy[x]-a-30
    right_gauss_bound = peaks.energy[x]+a
    data=data.loc[(data['energy'] >= left_gauss_bound) & (data['energy'] <= right_gauss_bound)]
    x_values_1 = np.asarray(data['energy'])
    y_values_1 = np.asarray(data['counts'])
    
    fit_function = input("Which function whould you like to use \n Gaussian (g), Gaussian+line (lg) or Gaussian+exponent (ag) > ")

    if fit_function=='g':
        p0 = [max(y_values_1),np.mean(x_values_1),np.std(x_values_1)]
        coeff, var_matrix = curve_fit(gaus, x_values_1, y_values_1, p0=p0)
        hist_gaus_fit = gaus(x_values_1, *coeff)

    if fit_function=='ag':
        tau=float(input('Any idea about tau of the peak? (roughtly equals to 1.44*half_life of the biggest peak) > '))
        p0 = [sum(y_values_1),np.mean(x_values_1),np.std(x_values_1),tau]
        coeff, var_matrix = curve_fit(alpha_gaus, x_values_1, y_values_1, p0=p0)
        hist_gaus_fit = alpha_gaus(x_values_1, *coeff)
        
    if fit_function=='lg':
        slope=(y_values_1[0]-y_values_1[len(y_values_1)-1])/(x_values_1[0]-x_values_1[len(y_values_1)-1])
        incr=y_values_1[0]-slope*x_values_1[0]
        p0 = [sum(y_values_1),np.mean(x_values_1),np.std(x_values_1),slope,incr]
        coeff, var_matrix = curve_fit(line_gaus, x_values_1, y_values_1, p0=p0)
        hist_gaus_fit = line_gaus(x_values_1, *coeff)
    
    fig2, ax2 = plt.subplots()
    plt.step(data.energy, data.counts,c='black')
    plt.plot(x_values_1,hist_gaus_fit, 'red', linewidth=2)
    ax2.tick_params(axis='both',direction='out', top=True,right=True)
    ax2.set_title('$^{48}Ca+^{242}Pu$', fontsize=20)
    ax2.set_ylabel('Counts', fontsize=18)
    ax2.set_xlabel('Energy, keV', fontsize=18)
    stderr_sigma_gaus = np.sqrt(var_matrix[2,2])
    stderr_mu_gaus = np.sqrt(var_matrix[1,1])
    a=coeff[1]
    b=stderr_mu_gaus
    c=coeff[2]
    d=stderr_sigma_gaus
    energy_and_counts=np.array(peaks.iloc[x])
    peak_parameters=str([energy_and_counts,a, b])
    plt.legend(('Data', 'Fit'), loc='upper left')
    title=file_name+' fit' +'\n$ \mu$='+str(round(a,3))+'$\pm$'+str(round(b,3))+'\n$ \sigma$='+str(round(c,3))+'$\pm$'+str(round(d,3))
    ax2.set_title(title)
    plt.ylim(min(data.counts), max(data.counts)+10)
    saving=input('Save the Fit? y/n > ')
    if saving =="y":
        plt.tight_layout()
        plt.savefig(file_name+"_Peak_"+str(x)+".pdf")
        file1 = open(file_name+"_Peak_"+str(x)+ ".txt","w")
        file1.write(peak_parameters)
        file1.close()
        plt.show()
    else:
        plt.show()
    cont = input("Repeat or fit another peak? y/n > ")
    if cont == "n":
        break
