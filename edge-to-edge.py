import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
	 
arf1gm130 = pd.read_excel (r'C:/Users/AG Bottanelli/Nextcloud2/Wong et al., 2021/ARF Paper 2021/Quantification/COPI - ARFs/Line Porfiles TOP View ARF1+COPI.xlsx')
'''
	Odd columns - intensity values 
    - I Ch1 (ARF): 4n+1 starting in 0 
    - I ch2 (marker): 4n-1 starting at 1 
  Even columns - distances 
'''
def distance_edge_to_edge(df):
    lineprofiles = np.array(df)
    lineprofiles = lineprofiles[2:,:]
    n = int(np.shape(lineprofiles)[1])
    numberprofiles = int(np.shape(lineprofiles)[1]/4)
    D=np.array(np.zeros(numberprofiles))  
    E=np.array(np.zeros(numberprofiles))
    #loop that extracts the distance and intensity values
    for i in range(0, n, 4): 
        j = int(abs(i/4)) 
        distance = lineprofiles[:,i] 
        ch1 = lineprofiles[:, i+1] 
        ch2 = lineprofiles[:,i+3]
        gradch1 = np.gradient(ch1, distance)
        gradch2 = np.gradient(ch2, distance) 
        np.nan_to_num(gradch1, 0)
        edgech11 =  gradch1.argmax() #Line profile drawn from ch2 to ch1  
        edgech12 = gradch1.argmin() #Line profile drawn from ch1 to ch2
        np.nan_to_num(gradch2, 0)
        edgech21 =  gradch2.argmax()      
        edgech22 = gradch2.argmin() 
        D[j] = (edgech21- edgech11)*30   
        E[j] = (edgech22- edgech12)*30
    #when drawing from channel channel 2 to channel 1, use edge 2 of channel 2 and edge 1 of channel 1. Vice versa
    return(D,E)
 
def plot_lineprofiles_grad(df): 
    lineprofiles = np.array(df)
    lineprofiles = lineprofiles[2:,:]
    n = int(np.shape(lineprofiles)[1])
    for i in range(0,n,4):
        j = int(abs(i/4)) 
        distance = lineprofiles[:,i]*1000
        ch1 = lineprofiles[:, i+1] 
        ch2 = lineprofiles[:,i+3]
        gradch1 = np.gradient(ch1, distance) 
        gradch1n = gradch1/np.nanmax(gradch1)        
        gradch2 = np.gradient(ch2, distance) 
        gradch2n = gradch2/np.nanmax(gradch2)
        ch1n = ch1/np.nanmax(ch1)
        ch2n = ch2/np.nanmax(ch2)
        plt.figure()
        plt.plot(distance, ch1n, label = 'GM signal')
        plt.plot(distance, ch2n, label = 'ARF1 signal')
        plt.plot(distance, gradch1n,'.', label = 'Gm gradient signal')      
        plt.plot(distance, gradch2n, '.', label = 'ARF1 gradient signal')
        plt.title('Line profile number: T={}'.format(j))
        plt.legend(loc="lower left", prop={'size': 7.5})
        plt.show()
 
Darf1gm130 = distance_edge_to_edge(arf1gm130)
plot_lineprofiles_grad(arf1gm130)
