# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:50:47 2020

@author: Ivan
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import peakutils
import re
import seaborn as sns
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#read the data and take only the meaningfull part

file_name=input('Please enter the isotope name to download the matrix>')
data = pd.read_csv('Matrix_'+file_name + '.dat',sep='	', decimal=',', float_precision='high', header=None)
data[data<0]=0

right_indexes=np.arange(len(data.index)-1,-1,-1)
right_indexes=right_indexes.astype(int)
data=data.set_index(right_indexes)
data=data.sort_index(axis=0)
x=[]

file_name_1=input('Please choose isotope_0 for colibration \n for example  Rn201_Peak_1 >')

with open(file_name_1 + '.txt', 'r') as f:
    first_line = f.readline()
    x=re.findall("(?<=[AZaz])?(?!\d*=)[0-9.+-]+",first_line)
peak_0=pd.DataFrame(x)

file_name_2=input('Please choose isotope_1 for colibration >')
with open(file_name_2 + '.txt', 'r') as f:
    first_line = f.readline()
    y=re.findall("(?<=[AZaz])?(?!\d*=)[0-9.+-]+",first_line)
peak_1=pd.DataFrame(y)

'''
while True:
    print('Please create the square M*M matrix for plotting, peak finding and colibration')
    x_min=float(input('x_min > '))
    x_max=float(input('x_max > '))
    y_min=float(input('y_min > '))
    y_max=float(input('y_max > '))
    x=np.arange(x_min,x_max,1)
    y=np.arange(y_min,y_max,1)
    data=data.take(x,axis=1)
    data=data.take(y,axis=0)
    cont = input("Repeat the plotting? y/n > ")
    fig1, ax1 = plt.subplots()
    ax1 = sns.heatmap(data, cmap="gist_gray_r", cbar_kws={'label': 'Counts'})
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0)
    for _, spine in ax1.spines.items():
        spine.set_visible(True)
    plt.show()
    if cont == "n":
           break

'''


a_list=[]
for i in range(0,len(data.index)):
    x_peaks=np.array(data.loc[i])
    threshold= 0.5
    minimum_distanse=10
    x_peaks= peakutils.indexes(x_peaks, thres=threshold, min_dist=minimum_distanse)
    #x_peaks= pd.DataFrame(data, index=x_peaks)
    a_list.append(x_peaks)
  
b_list=[]
for j in range(0,len(data.columns)):
    y_peaks=np.array(data.iloc[:,j])
    threshold= 0.5
    minimum_distanse=10
    y_peaks = peakutils.indexes(y_peaks, thres=threshold, min_dist=minimum_distanse)
    b_list.append(y_peaks)


df_x = pd.DataFrame(a_list)
df_y = pd.DataFrame(b_list)
df_x=df_x.reset_index()
df_y=df_y.reset_index()
df_x=df_x.melt(id_vars=["index"]).dropna(axis=0).sort_values(by=['index'])
df_y=df_y.melt(id_vars=["index"]).dropna(axis=0).sort_values(by=['index'])
df_x = df_x.rename(columns={'index': 'row', 'value': 'col'})
df_y = df_y.rename(columns={'index': 'col', 'value': 'row'})
df_x=df_x.reset_index()
df_y=df_y.reset_index()
df_x=df_x.drop(['index'], axis=1)
df_y=df_y.drop(['index'], axis=1)
df_x=df_x.drop(['variable'], axis=1)
df_y=df_y.drop(['variable'], axis=1)

columnsTitles=["row","col"]
df_y=df_y.reindex(columns=columnsTitles)

final_list=[]
for j in range(0,len(df_x)):
    a=np.array(df_x.loc[j])
    for i in range(0,len(df_y)):
        b=np.array(df_y.loc[i])
        equal_arrays=(a==b).all()
        if equal_arrays == True:
            final_list.append(b)

Final_peaks = pd.DataFrame(final_list)
Final_peaks.columns=["row", "col"]
row=np.array(Final_peaks.row)
col=np.array(Final_peaks.col)
counts=data.iloc[row, col]
#counts.columns = counts.columns.map(str)

counts=counts.groupby(counts.columns, axis=1).sum()
counts=counts.sort_index(axis=1)
counts = counts.groupby(level=0).last()

'''
while True:
    strip_1=int(input('Please inter the aproximate strip number of peak_0 referense > '))
    strip_2=int(input('Please inter the aproximate strip number of peak_1 referense > '))
    b=np.array(counts.columns)
    #equal_arrays_1=(strip_1==b).all()
    #equal_arrays_2=(strip_2==b).all()
    #print(strip_1)
    #print(strip_2)
    if strip_1 in b and strip_2 in b:
        y1=float(peak_0.loc[2])#counts.loc[:,strip_1].idxmax(axis=1)
        y2=float(peak_1.loc[2])#counts.loc[:,strip_2].idxmax(axis=1)
        x1=strip_1
        x2=strip_2
        a1=(y2-y1)/(x2-x1)
        b1=y1-a1*x1
        print(x1,y1)
        print(x2,y2)
        print(a1,b1)
        print("colibration finished successfully")
    else:
        if not strip_1 in b:
            print("No such string found for peak_0")
        if not strip_2 in b:
            print("No such string found for peak_1")
    cont = input("Repeat the colibration? y/n > ")
    if cont == "n":
           break
'''
y1=float(peak_0.loc[3])#counts.loc[:,strip_1].idxmax(axis=1)
y2=float(peak_1.loc[3])#counts.loc[:,strip_2].idxmax(axis=1)
a1=(y2-y1)/14
#first_element=y1-a1*98
element_0=y1-a1*53
element_1=element_0+a1*len(data.index)
colibration=np.arange(element_0,element_1,a1)
#colibration=a1*colibration
colibration=colibration.astype(int)
data=data.set_index(colibration)


while True:
    fig, ax = plt.subplots()
#data = data.pivot("Energy,Mev", "Strip number", "counts")
    ax = sns.heatmap(data, cmap="gist_gray_r", cbar_kws={'label': 'Counts'})
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    
    print('Strip range to plot the graph')
    x_min=float(input('x_min > '))
    x_max=float(input('x_max > '))
    
    print('Energy range to plot the graph (from 0 until lenght of the matrix)')
    y_min=float(input('y_min > '))
    y_max=float(input('y_max > '))
    plt.ylim(y_max, y_min)
    
    plt.xlim(x_min, x_max)
    
    plt.xlabel("Strip number", fontsize=18)
    plt.ylabel("Energy, kev",fontsize=18)
    ax.set_title("$^{48}Ca+^{242}Pu$", fontsize=20)
    while True:
        isotopes = input("add the name of the isotope on plot? y/n > ")
        if isotopes =="y":
           name=input("name of the isotope > ")
           x_pos= float(input("x position on the plot > "))
           y_pos=float(input("y position on the plot > "))
           plt.text(x_pos, y_pos, name)
        if isotopes == "n":
           break
    plt.tight_layout()
    plt.savefig("Heat_map_"+ file_name + ".pdf")
    plt.show()
    cont = input("Repeat the plotting? y/n > ")
    if cont == "n":
           break
#del x,y,file_name,file_name_1,file_name_2,first_line,threshold,minimum_distanse,y_peaks,x_peaks,a_list,b_list,Final_peaks,a,b,columnsTitles, df_x,df_y,equal_arrays, final_list, i,j,col,cont,isotopes,right_indexes,row,strip_1,strip_2,x1,x2,x_max,x_min,y1,y2,y_max,y_min
