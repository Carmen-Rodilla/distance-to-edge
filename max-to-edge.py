

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

arf1 = pd.read_excel (r'Absolute Pathway')
#INTRODUCE THE .xlsx files paths with the line profile information between the ‘’
#THE ORDER OF COLUMNS SHOULD BE DISTANCE CHANNEL ARF-INTENSITY CHANNEL ARF - D vesicle - I vesicle - again... 
 
'''
Odd columns - intensity values 
- arf: 4n+1 starting in 0 
- ves: 4n-1 starting at 1 
Even columns - distances 
'''


def distance_edge_to_max(df): 
    lineprofiles = np.array(df)
    lineprofiles = lineprofiles[2:,:]
    n = int(np.shape(lineprofiles)[1])
    numberprofiles = int(np.shape(lineprofiles)[1]/4)
    D=np.array(np.zeros(numberprofiles))
    #loop that extracts the distance and intensity values
    for i in range(0, n, 4): 
        j = int(abs(i/4)) 
        distance = lineprofiles[:,i] 
        arf = lineprofiles[:, i+1] 
        ves = lineprofiles[:,i+3]
        gradarf = np.gradient(arf, distance) 
        np.nan_to_num(gradarf, 0)
        edge =  gradarf.argmax() #Line profile drawn from vesicle to golgi       
        edge2 = gradarf.argmin() #Line profile drawn from golgi to vesicle
        np.nan_to_num(ves, 0)
        centerves = np.argmax(ves)
        D[j] = (edge - centerves)*30 #change edge to edge2 if line profile is drawn from golgi to vesicle 
    return(D)


def plot_lineprofiles(df):
    lineprofiles = np.array(df)
    lineprofiles = lineprofiles[2:,:]
    n = int(np.shape(lineprofiles)[1])
    for i in range(0,n,4):
        j = int(abs(i/4)) 
        distance = lineprofiles[:,i]*1000
        arf = lineprofiles[:, i+1] 
        ves = lineprofiles[:,i+3]
        gradarf = np.gradient(arf, distance) 
        gradarfn = gradarf/np.nanmax(gradarf)
        arfn = arf/np.nanmax(arf)
        vesn = ves/np.nanmax(ves)
        plt.figure()
        plt.plot(distance, arfn, label = 'ARF signal')
        plt.plot(distance, vesn, label = 'vesicle signal')
        plt.plot(distance, gradarfn, '.', label = 'Arf gradient signal')
        plt.title('Line profile number: T= {}'.format(j))
        plt.legend()
        plt.show()
        
        
        
Darf1 = distance_edge_to_max(arf1)
plot_lineprofiles(arf1)
