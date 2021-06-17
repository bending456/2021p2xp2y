import matplotlib.pylab as plt 
# routlines for analyzing odes
import sys
import pickle as pk
import numpy as np
import math
import analyzeGotran as ao
import subprocess as sb
from subprocess import PIPE
import shlex
import time
import ScriptRunner as SR
import matplotlib.pyplot as plt; plt.rcdefaults()
import scipy
import SArunner as SA
import matplotlib.cm as cm

## Contour Plot P2X vs. P2Y - TNFa mRNA, TNFa released, Degree of Chemotaxis 
## 1000 ATP 

change = np.array([0.1,0.5,1,2,3])
counter1 = 0
for i in np.arange(len(change)):
    counter2 = 0
    for j in np.arange(len(change)):
        data   = SR.gotranMicroglia(sim_time      = 300,
                                    ATP           = 1000,
                                    rhop2yc       = change[i],
                                    rhop2y12      = change[i],
                                    rhop2x4       = change[j],
                                    rhop2x7       = change[j],
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration29',
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter2 == 0:
            temppeakCa   = max(data[1])
            tempavgCa    = np.average(data[1][240:299])
            tempmRNATNF  = data[2][-1]-23.388667
            temppakt     = data[3][-1]
            tempdist     = data[4][-1]
            tempTNFae    = data[5][-1]-0.3387139
        else:
            temppeakCa   = np.append(temppeakCa,max(data[1]))
            tempavgCa    = np.append(tempavgCa,np.average(data[1][240:299]))
            tempmRNATNF  = np.append(tempmRNATNF,data[2][-1]-23.388667)
            temppakt     = np.append(temppakt,data[3][-1])
            tempdist     = np.append(tempdist,data[4][-1])
            tempTNFae    = np.append(tempTNFae,data[5][-1]-0.3387139)

        counter2 = counter2 + 1 
    
    if counter1 == 0:
        Capeak = temppeakCa
        Caavg  = tempavgCa
        mRNA   = tempmRNATNF
        pAkt   = temppakt
        dist   = tempdist
        TNFa   = tempTNFae
    else:
        Capeak = np.vstack([Capeak,temppeakCa])
        Caavg  = np.vstack([Caavg,tempavgCa])
        mRNA   = np.vstack([mRNA,tempmRNATNF])
        pAkt   = np.vstack([pAkt,temppakt])
        dist   = np.vstack([dist,tempdist])
        TNFa   = np.vstack([TNFa,tempTNFae])
        
    counter1 = counter1 + 1

standard = SR.gotranMicroglia(sim_time      = 300,
                              ATP           = 10,
                              rhop2yc       = 1,
                              rhop2x4       = 1,
                              rhop2x7       = 1,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration29',
                              data_name2    = 'mRNA_TNF',
                              data_name3    = 'pAkt',
                              data_name4    = 'S',
                              data_name5    = 'TNFae',
                              removePickle  = 1,
                              timePrint     = 0)
Caapex = max(standard[1]-100.038)
Catavg = np.average(standard[1][240:299])-100.038
mRNAstd = standard[2][-1]-23.388667
paktstd = standard[3][-1]
diststd = standard[4][-1]
TNFaestd = standard[5][-1]-0.3387139  
    
data13 = (Capeak) #(Capeak-Capeak[3][3])/Capeak[6][6]
data14 = (Caavg) #(Caavg-Caavg[3][3])/Caavg[6][6]
data15 = (mRNA) #(mRNA-mRNA[3][3])/mRNA[6][6]
data16 = (pAkt) #(pAkt-pAkt[3][3])/pAkt[6][6]
data17 = (dist) #(dist-dist[3][3])/dist[6][6]
data18 = (TNFa) #(TNFa-TNFa[3][3])/TNFa[6][6]

data171 = data17/data17[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data171, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5) #(distmin,distmax,diststep)
CS = plt.contour(data171, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Degree of Chemotaxis (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Expression (in log scale)",fontsize=15)
plt.ylabel("P2Y Expression (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YDist1000ATP.png')

data181 = data18/data18[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data181, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5) #(tnfamin,tnfamax,tnfastep)
CS = plt.contour(data181, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Released TNFa (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Expression (in log scale)",fontsize=15)
plt.ylabel("P2Y Expression (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YTNFa1000ATP.png')

## 10 ATP 

change = np.array([0.1,0.5,1,2,3])
counter1 = 0
for i in np.arange(len(change)):
    counter2 = 0
    for j in np.arange(len(change)):
        data   = SR.gotranMicroglia(sim_time      = 300,
                                    ATP           = 10,
                                    rhop2yc       = change[i],
                                    rhop2y12      = change[i],
                                    rhop2x4       = change[j],
                                    rhop2x7       = change[j],
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration29',
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter2 == 0:
            temppeakCa   = max(data[1])
            tempavgCa    = np.average(data[1][240:299])
            tempmRNATNF  = data[2][-1]-23.388667
            temppakt     = data[3][-1]
            tempdist     = data[4][-1]
            tempTNFae    = data[5][-1]-0.3387139
        else:
            temppeakCa   = np.append(temppeakCa,max(data[1]))
            tempavgCa    = np.append(tempavgCa,np.average(data[1][240:299]))
            tempmRNATNF  = np.append(tempmRNATNF,data[2][-1]-23.388667)
            temppakt     = np.append(temppakt,data[3][-1])
            tempdist     = np.append(tempdist,data[4][-1])
            tempTNFae    = np.append(tempTNFae,data[5][-1]-0.3387139)

        counter2 = counter2 + 1 
    
    if counter1 == 0:
        Capeak = temppeakCa
        Caavg  = tempavgCa
        mRNA   = tempmRNATNF
        pAkt   = temppakt
        dist   = tempdist
        TNFa   = tempTNFae
    else:
        Capeak = np.vstack([Capeak,temppeakCa])
        Caavg  = np.vstack([Caavg,tempavgCa])
        mRNA   = np.vstack([mRNA,tempmRNATNF])
        pAkt   = np.vstack([pAkt,temppakt])
        dist   = np.vstack([dist,tempdist])
        TNFa   = np.vstack([TNFa,tempTNFae])
        
    counter1 = counter1 + 1

standard = SR.gotranMicroglia(sim_time      = 300,
                              ATP           = 10,
                              rhop2yc       = 1,
                              rhop2x4       = 1,
                              rhop2x7       = 1,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration29',
                              data_name2    = 'mRNA_TNF',
                              data_name3    = 'pAkt',
                              data_name4    = 'S',
                              data_name5    = 'TNFae',
                              removePickle  = 1,
                              timePrint     = 0)
Caapex = max(standard[1]-100.038)
Catavg = np.average(standard[1][240:299])-100.038
mRNAstd = standard[2][-1]-23.388667
paktstd = standard[3][-1]
diststd = standard[4][-1]
TNFaestd = standard[5][-1]-0.3387139

data1 = (Capeak)
data2 = (Caavg)
data3 = (mRNA)
data4 = (pAkt)
data5 = (dist)
data6 = (TNFa)

data51 = data5/data5[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data51, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5)
CS = plt.contour(data51, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Degree of Chemotaxis (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Expression (in log scale)",fontsize=15)
plt.ylabel("P2Y Expression (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YDist10ATP.png')

data61 = data6/data6[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data61, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5)
CS = plt.contour(data61, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Released TNFa (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Expression (in log scale)",fontsize=15)
plt.ylabel("P2Y Expression (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YTNFa10ATP.png')

## 100 ATP 
change = np.array([0.1,0.5,1,2,3])
counter1 = 0
for i in np.arange(len(change)):
    counter2 = 0
    for j in np.arange(len(change)):
        data   = SR.gotranMicroglia(sim_time      = 300,
                                    ATP           = 100,
                                    rhop2yc       = change[i],
                                    rhop2x4       = change[j],
                                    rhop2x7       = change[j],
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration27',
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter2 == 0:
            temppeakCa   = max(data[1])
            tempavgCa    = np.average(data[1][240:299])
            tempmRNATNF  = data[2][-1]-23.388667
            temppakt     = data[3][-1]
            tempdist     = data[4][-1]
            tempTNFae    = data[5][-1]-0.3387139
        else:
            temppeakCa   = np.append(temppeakCa,max(data[1]))
            tempavgCa    = np.append(tempavgCa,np.average(data[1][240:299]))
            tempmRNATNF  = np.append(tempmRNATNF,data[2][-1]-23.388667)
            temppakt     = np.append(temppakt,data[3][-1])
            tempdist     = np.append(tempdist,data[4][-1])
            tempTNFae    = np.append(tempTNFae,data[5][-1]-0.3387139)

        counter2 = counter2 + 1 
    
    if counter1 == 0:
        Capeak = temppeakCa
        Caavg  = tempavgCa
        mRNA   = tempmRNATNF
        pAkt   = temppakt
        dist   = tempdist
        TNFa   = tempTNFae
    else:
        Capeak = np.vstack([Capeak,temppeakCa])
        Caavg  = np.vstack([Caavg,tempavgCa])
        mRNA   = np.vstack([mRNA,tempmRNATNF])
        pAkt   = np.vstack([pAkt,temppakt])
        dist   = np.vstack([dist,tempdist])
        TNFa   = np.vstack([TNFa,tempTNFae])
        
    counter1 = counter1 + 1
    
data7 = (Capeak)
data8 = (Caavg) 
data9 = (mRNA) 
data10 = (pAkt) 
data11 = (dist) 
data12 = (TNFa) 

data111 = data11/data11[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data111, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5) #(distmin,distmax,diststep)
CS = plt.contour(data111, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Degree of Chemotaxis (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Density (in log scale)",fontsize=15)
plt.ylabel("P2Y Density (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YDist10ATP.png')

data121 = data12/data12[4][4]
plt.figure(figsize=(7,8),dpi=150)
plt.tick_params(labelsize=12)
im = plt.imshow(data121, interpolation='bilinear', origin='lower',
                cmap=cm.summer, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)),aspect='auto')
levels = scipy.linspace(0,1,5) #(tnfamin,tnfamax,tnfastep)
CS = plt.contour(data121, levels, origin='lower',colors='k', 
                linewidths=2, extent=(np.log(0.1)/np.log(10), np.log(3)/np.log(10), np.log(0.1)/np.log(10), np.log(3)/np.log(10)))

CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
CBI.set_label('Released TNFa (Scale to Max)', fontsize=12)
plt.clabel(CS, inline=1, fontsize=12)
plt.xlabel("P2X Density (in log scale)",fontsize=15)
plt.ylabel("P2Y Density (in log scale)",fontsize=15)
plt.tight_layout()
plt.savefig('ContP2XP2YTNFa10ATP.png')


import datetime
import numpy as np
import scipy as sp
import scipy.fftpack
import pandas as pd

from spectrum import *
from scipy.interpolate import spline

## Ca transients 

ATP = np.array([1000])
time = 600
counter = 0
p2x4 = np.array([1])
p2x7 = np.array([0,1])
p2yc = np.array([1])
for i in np.arange(len(ATP)): # ATP
    for j in np.arange(len(p2x4)):
        for k in np.arange(len(p2x7)):
            for l in np.arange(len(p2yc)):
                if p2x4[j] + p2x7[k] + p2yc[l] == 0:
                    continue
                else:
                    data   = SR.gotranMicroglia(sim_time      = time,
                                                ATP           = ATP[i],
                                                output_name   = 'test1',
                                                ode_file_name = 'p2xp2yMigration29',
                                                rhop2x4       = p2x4[j],
                                                rhop2x7       = p2x7[k],
                                                rhop2yc       = p2yc[l],
                                                data_name2    = 'TNFae',
                                                data_name3    = 'TNFac',
                                                data_name4    = 'S',
                                                data_name5    = 'NFATNn',
                                                data_name6    = 'PI3Ka',
                                                data_name7    = 'pAkt',
                                                data_name8    = 'Ca4_CN',
                                                data_name9    = 'CaMCN',
                                                EleSwitch     = 0,
                                                DegSwitch     = 1,
                                                removePickle  = 1)
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        te      = data[2]
                        tc      = data[3]
                        dist    = data[4]
                        nf      = data[5]
                        p3      = data[6]
                        pa      = data[7]
                        cn      = 0.1*data[8] + data[9]
                        atp     = ATP[i]
                        p4      = p2x4[j]
                        p7      = p2x7[k]
                        py      = p2yc[l]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        te      = np.vstack([te,data[2]])
                        tc      = np.vstack([tc,data[3]])
                        dist    = np.vstack([dist,data[4]])
                        nf      = np.vstack([nf,data[5]])
                        p3      = np.vstack([p3,data[6]])
                        pa      = np.vstack([pa,data[7]])
                        cn      = np.vstack([cn,data[8]])
                        atp     = np.append(atp,ATP[i])
                        p4      = np.append(p4,p2x4[j])
                        p7      = np.append(p7,p2x7[k])
                        py      = np.append(py,p2yc[l])
                              
                    counter = counter + 1
                    
                    
noP2X7 = np.array([100.00,104.00,168.80,128.53,121.60,118.93,123.20,116.27,108.80,109.87,107.73,110.67,
                   105.33,109.07,102.93,101.33,105.60,104.00,106.13,104.00,113.87,105.33,104.53,103.20,
                   109.33,106.67,103.20,102.67,104.27,102.40]) 
time1 =  np.array([0.00,0.42,0.61,0.98,1.10,1.33,1.44,1.86,2.05,2.08,2.35,2.54,2.73,2.92,3.33,3.64,3.79,3.98,4.09,4.17,
                   4.58,4.89,5.15,5.42,5.68,5.95,6.14,6.52,6.63,6.86])

lowATP = np.array([100.00,100.53,136.80,178.93,174.93,130.13,107.47,100.53,100.27,100.80,100.27,
                   101.60,101.87,100.53,100.80,100.80,101.33,100.53,102.93])
time2  = np.array([0.00,0.38,0.49,1.02,1.21,1.67,2.16,2.77,3.11,4.24,4.32,4.66,4.89,5.08,5.30,5.49,5.91,6.10,6.40])

All    = np.array([100.00,100.00,175.73,166.13,156.27,161.60,153.60,162.93,155.47,164.53,162.40,163.47,154.40,160.27,
                   151.73,161.60,136.27,157.33,137.07,131.73,134.13,134.13,157.87,131.73,146.93,128.80,129.60,129.07,
                   131.47,151.73,125.87,133.33,125.87,149.60,126.40])
time3  = np.array([0.00,0.42,0.61,0.72,1.06,1.17,1.33,1.48,1.93,2.01,
                   2.16,2.20,2.42,2.50,2.69,2.80,3.18,3.41,3.60,3.83,
                   3.94,4.09,4.28,4.70,5.27,5.80,5.91,6.10,6.29,6.52,
                   7.08,7.50,7.92,8.26,8.79])

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
plt.plot((time3-0.42)*60,All,'r--',alpha=0.65,lw=2,label='Expt. Control')
plt.plot((time1-0.42)*60,noP2X7,'b--',alpha=0.65,lw=2,label='Expt. P2X7 inhibited') 
plt.plot(dura,Ca[1],'r',lw=1.5,alpha=1,label='Control')
plt.plot(dura,Ca[0],'b-',lw=1.5,label='P2X7 KO')
plt.plot([5,595],[97,97],'k',lw=5,label='ATP')
plt.xlim([-20,520])

plt.xlabel('Time (sec)',fontsize=12)
plt.ylabel('$[Ca^{2+}]_i$ (nM)')
plt.legend(fontsize=10,loc=0)
plt.tight_layout()
plt.savefig('CaP2XvsP2Y.png')


ATP = np.array([100,1000])
time = 600
counter = 0
p2x4 = np.array([0])
p2x7 = np.array([0,1])
p2yc = np.array([1])
for i in np.arange(len(ATP)): # ATP
    for j in np.arange(len(p2x4)):
        for k in np.arange(len(p2x7)):
            for l in np.arange(len(p2yc)):
                if p2x4[j] + p2x7[k] + p2yc[l] == 0:
                    continue
                else:
                    data   = SR.gotranMicroglia(sim_time      = time,
                                                ATP           = ATP[i],
                                                output_name   = 'test1',
                                                ode_file_name = 'p2xp2yMigration29',
                                                rhop2x4       = 0,
                                                rhop2x7       = 0,
                                                rhop2yc       = p2yc[l],
                                                data_name2    = 'TNFae',
                                                data_name3    = 'TNFac',
                                                data_name4    = 'S',
                                                data_name5    = 'NFATNn',
                                                data_name6    = 'PI3Ka',
                                                data_name7    = 'pAkt',
                                                data_name8    = 'Ca4_CN',
                                                data_name9    = 'CaMCN',
                                                EleSwitch     = 0,
                                                DegSwitch     = p2x7[k],
                                                removePickle  = 1)
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        te      = data[2]
                        tc      = data[3]
                        dist    = data[4]
                        nf      = data[5]
                        p3      = data[6]
                        pa      = data[7]
                        cn      = 0.1*data[8] + data[9]
                        atp     = ATP[i]
                        p4      = p2x4[j]
                        p7      = p2x7[k]
                        py      = p2yc[l]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        te      = np.vstack([te,data[2]])
                        tc      = np.vstack([tc,data[3]])
                        dist    = np.vstack([dist,data[4]])
                        nf      = np.vstack([nf,data[5]])
                        p3      = np.vstack([p3,data[6]])
                        pa      = np.vstack([pa,data[7]])
                        cn      = np.vstack([cn,data[8]])
                        atp     = np.append(atp,ATP[i])
                        p4      = np.append(p4,p2x4[j])
                        p7      = np.append(p7,p2x7[k])
                        py      = np.append(py,p2yc[l])
                              
                    counter = counter + 1
                    
                    
                    
plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
plt.plot(dura,Ca[0],'k--',lw=1.5,alpha=0.75,label='100 uM')
plt.plot(dura,Ca[2],'k-',lw=1.5,alpha=1,label='1000 uM')
plt.plot([5,595],[97,97],'k',lw=5,label='ATP')
plt.xlim([-20,520])
plt.ylim([90,160])

plt.xlabel('Time (sec)',fontsize=12)
plt.ylabel('$[Ca^{2+}]_i$ (nM)')
plt.legend(fontsize=10,loc=0)
plt.tight_layout()
plt.savefig('CaP2YonlynoDeg.png')

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
plt.plot(dura,Ca[1],'k--',lw=1.5,alpha=0.75,label='100 uM')
plt.plot(dura,Ca[3],'k-',lw=1.5,alpha=1,label='1000 uM')
plt.plot([5,595],[97,97],'k',lw=5,label='ATP')
plt.xlim([-20,520])
plt.ylim([90,160])


plt.xlabel('Time (sec)',fontsize=12)
plt.ylabel('$[Ca^{2+}]_i$ (nM)')
plt.legend(fontsize=10,loc=0)
plt.tight_layout()
plt.savefig('CaP2YonlywDeg.png')


ATP = np.array([10,100,500,1000])
time = 600
counter = 0
p2x4 = np.array([1])
p2x7 = np.array([1])
p2yc = np.array([1])
for i in np.arange(len(ATP)): # ATP
    for j in np.arange(len(p2x4)):
        for k in np.arange(len(p2x7)):
            for l in np.arange(len(p2yc)):
                if p2x4[j] + p2x7[k] + p2yc[l] == 0:
                    continue
                else:
                    data   = SR.gotranMicroglia(sim_time      = time,
                                                ATP           = ATP[i],
                                                output_name   = 'test1',
                                                ode_file_name = 'p2xp2yMigration28',
                                                rhop2x4       = 1,
                                                rhop2x7       = 1,
                                                rhop2yc       = 1,
                                                EleSwitch     = 0,
                                                DegSwitch     = 1,
                                                removePickle  = 1)
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        atp     = ATP[i]
                        p4      = p2x4[j]
                        p7      = p2x7[k]
                        py      = p2yc[l]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        atp     = np.append(atp,ATP[i])
                        p4      = np.append(p4,p2x4[j])
                        p7      = np.append(p7,p2x7[k])
                        py      = np.append(py,p2yc[l])
                              
                    counter = counter + 1
                    

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
plt.plot(dura,Ca[0],'k--',lw=1.5,alpha=0.5,label='10 uM')
plt.plot(dura,Ca[1],'k-.',lw=1.5,alpha=0.65,label='100 uM')
plt.plot(dura,Ca[2],'k:',lw=1.5,alpha=1,label='500 uM')
plt.plot(dura,Ca[3],'k-',lw=1.5,alpha=1,label='1000 uM')

plt.plot([5,595],[97,97],'k',lw=5,label='ATP')
plt.xlim([-20,520])

plt.xlabel('Time (sec)',fontsize=12)
plt.ylabel('$[Ca^{2+}]_i$ (nM)')
plt.legend(fontsize=10,loc=0)
plt.tight_layout()
plt.savefig('CavsATP.png')


## Gprot vs. DAG

data10   = SR.gotranMicroglia(sim_time      = 10000,
                              ATP           = 10,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration29',
                              data_name2    = 'GaGTP',
                              data_name3    = 'DAG',
                              DegSwitch     = 0,
                              removePickle  = 1)

data100   = SR.gotranMicroglia(sim_time      = 10000,
                              ATP           = 100,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration29',
                              data_name2    = 'GaGTP',
                              data_name3    = 'DAG',
                              DegSwitch     = 0,
                              removePickle  = 1)

data1000   = SR.gotranMicroglia(sim_time      = 10000,
                              ATP           = 1000,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration29',
                              data_name2    = 'GaGTP',
                              data_name3    = 'DAG',
                              DegSwitch     = 0,
                              removePickle  = 1)

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
#plt.plot(data10[2],data10[3],'bo',alpha=1,label='10 uM')
plt.plot(data100[2],data100[3],'rs',alpha=0.5,label='100 uM')
plt.plot(data1000[2],data1000[3],'gx',alpha=0.5,label='1000 uM')
plt.xlabel("Active G protein [a.u.]",fontsize=12)
plt.ylabel("IP3 & DAG [a.u.]",fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()
plt.savefig('withoutDeg.png')

data100n   = SR.gotranMicroglia(sim_time      = 10000,
                              ATP           = 100,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration27',
                              #rhop2y12      = p2yc[l],
                              data_name2    = 'GaGTP',
                              data_name3    = 'DAG',
                              DegSwitch     = 1,
                              removePickle  = 1)

data1000n   = SR.gotranMicroglia(sim_time      = 10000,
                              ATP           = 1000,
                              output_name   = 'test1',
                              ode_file_name = 'p2xp2yMigration27',
                              #rhop2y12      = p2yc[l],
                              data_name2    = 'GaGTP',
                              data_name3    = 'DAG',
                              DegSwitch     = 1,
                              removePickle  = 1)

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
#plt.plot(data10[2],data10[3],'bo',alpha=1,label='10 uM')
plt.plot(data100n[2],data100n[3],'rs',alpha=0.5,label='100 uM')
plt.plot(data1000n[2],data1000n[3],'gx',alpha=0.5,label='1000 uM')
plt.xlabel("Active G protein [a.u.]",fontsize=12)
plt.ylabel("IP3 & DAG [a.u.]",fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()
plt.savefig('withDeg.png')


p2x4 = np.array([0,1,9])
p2x7 = np.array([0,1])
p2y1 = 1
scale = 10000
###### No ATP ####################################
dtime     = SR.gotranMicroglia(sim_time      = 10800,
                               ATP           = 3000,
                               output_name   = 'test1',
                               ode_file_name = 'p2xp2yMigration29',
                               data_name2    = 'DNA_TNF',
                               data_name3    = 'DNATNF',
                               data_name4    = 'mRNA_TNF',
                               data_name5    = 'TNFae',
                               DegSwitch     = 0,
                               removePickle  = 1)

d1     = SR.gotranMicroglia(sim_time      = 10800,
                            ATP           = 1000,
                            output_name   = 'test1',
                               ode_file_name = 'p2xp2yMigration29',
                               data_name2    = 'DNA_TNF',
                               data_name3    = 'DNATNF',
                               data_name4    = 'mRNA_TNF',
                               data_name5    = 'TNFae',
                               DegSwitch     = 0,
                               removePickle  = 1)

d02     = SR.gotranMicroglia(sim_time      = 10800,
                             ATP           = 300,
                             output_name   = 'test1',
                               ode_file_name = 'p2xp2yMigration29',
                               data_name2    = 'DNA_TNF',
                               data_name3    = 'DNATNF',
                               data_name4    = 'mRNA_TNF',
                               data_name5    = 'TNFae',
                               DegSwitch     = 0,
                               removePickle  = 1)

d01     = SR.gotranMicroglia(sim_time      = 10800,
                            ATP           = 100,
                            output_name   = 'test1',
                               ode_file_name = 'p2xp2yMigration29',
                               data_name2    = 'DNA_TNF',
                               data_name3    = 'DNATNF',
                               data_name4    = 'mRNA_TNF',
                               data_name5    = 'TNFae',
                               DegSwitch     = 0,
                               removePickle  = 1)

d001     = SR.gotranMicroglia(sim_time      = 10800,
                            ATP           = 10,
                            output_name   = 'test1',
                               ode_file_name = 'p2xp2yMigration29',
                               data_name2    = 'DNA_TNF',
                               data_name3    = 'DNATNF',
                               data_name4    = 'mRNA_TNF',
                               data_name5    = 'TNFae',
                               DegSwitch     = 0,
                               removePickle  = 1)

hrs = np.array([0,1,2,3])
Hide = np.array([0,0.054,0.54,1])

calctnfamrna0 = np.array([dtime[5][0],dtime[5][3599],dtime[5][7199],dtime[5][10799]])-dtime[5][0]
calctnfamrna1 = (calctnfamrna0)/max(calctnfamrna0)

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12)
plt.plot(hrs,calctnfamrna1,'b-s',label="Model- 3 mM ATP")
plt.plot(hrs,Hide,'r--s',label="Hide et al. - 3 mM ATP")

plt.xlim(-0.2,3.2)
plt.ylim(-0.05,1.05)

plt.legend(loc=0,fontsize=10)
plt.xlabel("Time [hr]",fontsize=12)
plt.ylabel("TNFa [fraction of max]",fontsize=12)
plt.tight_layout()
plt.savefig('tnftime.png')

calctnfamrna2 = np.array([d001[5][-1],d01[5][-1],d02[5][-1],d1[5][-1]])
calctnfamran3 = (calctnfamrna2-min(calctnfamrna2))/max(calctnfamrna2-min(calctnfamrna2))

ATP1 = np.array([0.01,0.1,0.3,1])
ATP2 = np.array([0.01,0.1,0.3,1])
Hide2 = np.array([0.04,0.26,0.35,1])

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12)
plt.semilogx(ATP2,calctnfamran3,'b-s',label="Model - 3 hr")
plt.semilogx(ATP1,Hide2,'r--s',label="Hide et al. - 3 hr")

plt.xlim(0.009,1.2)
#plt.ylim(0,1.05)

plt.legend(loc=0,fontsize=10)
plt.xlabel("ATP [mM]",fontsize=12)
plt.ylabel("TNFa [fraction of max]",fontsize=12)
plt.tight_layout()
plt.savefig('tnfdoselog.png')


plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12)
plt.plot(ATP2,calctnfamran3,'b-s',label="Model - 3 hr")
plt.plot(ATP1,Hide2,'r--s',label="Hide et al. - 3 hr")

plt.legend(loc=0,fontsize=10)
plt.xlabel("ATP (mM)",fontsize=12)
plt.ylabel("Synthesized TNFa (fraction of max)",fontsize=12)
plt.tight_layout()
plt.savefig('tnfdose.png')

## NFAT validation 

###### No ATP ####################################
d0     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 3000,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)


d00     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 0,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)

d100     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 100,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)

d500     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 500,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)
d1000     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 1000,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)

d3000     = SR.gotranMicroglia(sim_time      = 3600,
                            ATP           = 3000,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration29',
                            data_name2    = 'NFATpc',
                            data_name3    = 'NFATpn',
                            data_name4    = 'NFATNn',
                            data_name5    = 'Ca4_CN',
                            DegSwitch     = 0,
                            removePickle  = 1)

LitNFATr = np.array([0,0.26,0.95,1,0.95]) # in fraction
Littime = np.array([0,1,15,30,60]) # in min
CalcData = np.array([d0[4][0],d0[4][59],d0[4][899],d0[4][1799],d0[4][3599]])-d0[4][0]
ProcData = CalcData/max(CalcData)

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12)
plt.plot(Littime,ProcData,'b-s',label="Model at 3 mM [ATP]")
plt.plot(Littime,LitNFATr,'r--s',label="Ferrari et al., at 3 mM [ATP]")

plt.xlim(-2,62)
plt.ylim(-0.05,1.05)

plt.legend(loc=0,fontsize=10)
plt.xlabel("Time [min]",fontsize=12)
plt.ylabel("Dephosphorylated NFAT [Scale]",fontsize=12)
plt.tight_layout()
plt.savefig('nfattime.png')

ATPmM1 = np.array([0,0.1,0.5,1,3])
ATPmM2 = np.array([0,0.1,0.5,1,3,5]) 
LitNFATr = np.array([0,0.05,0.1,0.5,1,0.3])
CalcData2 = np.array([d00[4][-1],d100[4][-1],d500[4][-1],d1000[4][-1],d3000[4][-1]])-d00[4][-1]
ProcData2 = CalcData2/max(CalcData2)

plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12)
plt.plot(ATPmM1,ProcData2,'b-s',label="Model for 15 min")
plt.plot(ATPmM2,LitNFATr,'r--s',label="Ferrari et. al., for 15 min")

plt.xlim(-0.1,5.1)
plt.ylim(-0.05,1.05)

plt.legend(loc=0,fontsize=10)
plt.xlabel("ATP [mM]",fontsize=12)
plt.ylabel("Dephosphorylated NFAT [Scale]",fontsize=12)
plt.tight_layout()
plt.savefig('nfatdose.png')



## Akt validation 

ATP = np.array([10,50,100,1000])
time = 600
counter = 0
p2x4 = np.array([1])
p2x7 = np.array([1])
p2yc = np.array([1])
for i in np.arange(len(ATP)): # ATP
    for j in np.arange(len(p2x4)):
        for k in np.arange(len(p2x7)):
            for l in np.arange(len(p2yc)):
                if p2x4[j] + p2x7[k] + p2yc[l] == 0:
                    continue
                elif p2x4[j] == 0 and p2x7[k] == 0:
                    continue
                elif p2x4[j] == 0 and p2yc[l] == 0:
                    continue
                elif p2x7[k] == 0 and p2yc[l] == 0:
                    continue
                else: 
                    data   = SR.gotranMicroglia(sim_time      = time,
                                        ATP           = ATP[i],
                                        output_name   = 'test1',
                                        ode_file_name = 'p2xp2yMigration29',
                                        rhop2x4       = p2x4[j],
                                        rhop2x7       = p2x7[k],
                                        rhop2yc       = p2yc[l],
                                        data_name2    = 'pAkt',
                                        removePickle  = 1)
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        pAkt    = data[2]
                        atp     = ATP[i]
                        p4      = p2x4[j]
                        p7      = p2x7[k]
                        py      = p2yc[l]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        pAkt    = np.vstack([pAkt,data[2]])
                        atp     = np.append(atp,ATP[i])
                        p4      = np.append(p4,p2x4[j])
                        p7      = np.append(p7,p2x7[k])
                        py      = np.append(py,p2yc[l])
                              
                    counter = counter + 1
                    
                    
plt.figure(figsize=(6,4),dpi=150)
plt.tick_params(labelsize=12,direction='in')
#plt.title('A',fontsize=12,fontweight='bold',loc='left')
time = [0,1,2,5,10]
model = np.array([pAkt[1][0],pAkt[1][59],pAkt[1][119],pAkt[1][299],pAkt[1][599]])/pAkt[1][599]*100
#plt.plot(time,model,'b-s',linewidth=1.5,alpha=1,label="Model")
plt.plot(time,[0,30,70,95,100],'r--s',linewidth=1.5,alpha=1,label='Ohsawa $et$ $al.$ - 50 uM')
plt.tight_layout()

#plt.title('B',fontsize=12,fontweight='bold',loc='left')
plt.plot(dura/60,pAkt[0]/pAkt[1][599]*100,'k:',linewidth=2.5,alpha=0.25,label="10 uM")
plt.plot(dura/60,pAkt[1]/pAkt[1][599]*100,'k-.',linewidth=2.5,alpha=0.5,label="50 uM")
plt.plot(dura/60,pAkt[2]/pAkt[1][599]*100,'k--',linewidth=2.5,alpha=0.75,label="100 uM")
plt.plot(dura/60,pAkt[3]/pAkt[1][599]*100,'k-',linewidth=2.5,alpha=1,label="1000 uM")
plt.legend(loc=0)
plt.xlabel("Time (min)",fontsize=12)
plt.ylabel("Akt Phosphorylation (relative %)",fontsize=12)
#plt.xlim([-0.1,6.2])
plt.tight_layout()
plt.subplots_adjust(wspace=0.25)
plt.savefig('pAkt.png')

ATP = np.array([10,50,100,1000])
time = 300
counter = 0
p2x4 = np.array([0,1])
p2x7 = np.array([0,1])
p2yc = np.array([0,1])
for i in np.arange(len(ATP)): # ATP
    for j in np.arange(len(p2x4)):
        for k in np.arange(len(p2x7)):
            for l in np.arange(len(p2yc)):
                if p2x4[j] + p2x7[k] + p2yc[l] == 0:
                    continue
                elif p2x4[j] == 0 and p2x7[k] == 0:
                    continue
                elif p2x4[j] == 0 and p2yc[l] == 0:
                    continue
                elif p2x7[k] == 0 and p2yc[l] == 0:
                    continue
                else: 
                    data   = SR.gotranMicroglia(sim_time      = time,
                                            ATP           = ATP[i],
                                            output_name   = 'test1',
                                            ode_file_name = 'p2xp2yMigration29',
                                            rhop2x4       = p2x4[j],
                                            rhop2x7       = p2x7[k],
                                            rhop2y12      = p2yc[l],
                                            data_name2    = 'Q12_ptxf',
                                            data_name3    = 'Gb',
                                            data_name4    = 'PI3Ka',
                                            data_name5    = 'pAkt',
                                            data_name6    = 'p2y12a',
                                            data_name7    = 'S',
                                            DegSwitch     = 1,
                                            removePickle  = 1)
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        Q12     = data[2]
                        Gb      = data[3]
                        pi3ka   = data[4]
                        pakt    = data[5]
                        p2y12   = data[6]
                        S2      = data[7]
                        atp     = ATP[i]
                        p4      = p2x4[j]
                        p7      = p2x7[k]
                        p12     = p2yc[l]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        Q12     = np.vstack([Q12,data[2]])
                        Gb      = np.vstack([Gb,data[3]])
                        pi3ka   = np.vstack([pi3ka,data[4]])
                        pakt    = np.vstack([pakt,data[5]])
                        p2y12   = np.vstack([p2y12,data[6]])
                        S2      = np.vstack([S2,data[7]])
                        atp     = np.append(atp,ATP[i])
                        p4      = np.append(p4,p2x4[j])
                        p7      = np.append(p7,p2x7[k])
                        p12     = np.append(p12,p2yc[l])
                              
                    counter = counter + 1
    
data = np.array([pakt[7][300],pakt[5][300],pakt[4][300],pakt[6][300]])
ratio = data/max(data)
expratio = np.array([100,90,43,8])/100
objects = ('Control','P2X7 -/-','P2X4 -/-','P2Y -/-')

y_pos = np.arange(len(objects))


bar_width = 0.3
plt.figure(figsize=(12,4),dpi=150)
ax = plt.subplot(1,2,1)
plt.title('A',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-0.15, expratio, bar_width, align='center', color='r', alpha=1, edgecolor='black', label='Exp.')
plt.bar(y_pos+0.15, ratio, bar_width, align='center', color='b', alpha=1, edgecolor='black', label='Model')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('Akt phosphorylation [relative scale]',fontsize=12)
plt.legend(loc=0,fontsize=12)
plt.tight_layout()

data2 = np.array([S2[7][300],S2[5][300],S2[4][300],S2[6][300]])
ratio = data2/max(data2)
expratio = np.array([100,70,46,19])/100
objects = ('Control','P2X7 -/-','P2X4 -/-','P2Y -/-')

y_pos = np.arange(len(objects))
ax = plt.subplot(1,2,2)
plt.title('B',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-0.15, expratio, bar_width, align='center', color='r', alpha=1, edgecolor='black', label='Exp.')
plt.bar(y_pos+0.15, ratio, bar_width, align='center', color='b', alpha=1, edgecolor='black', label='Model')
plt.tick_params(labelsize=12,direction='in')
#plt.ylim([0.5,2.85])
plt.xticks(y_pos, objects)
plt.ylabel('Mig. Dist./Chemotaxis [relative scale]',fontsize=12)
plt.legend(loc=0,fontsize=12)
plt.tight_layout()

plt.subplots_adjust(wspace=0.25)
plt.savefig("pAktandMigration.png")


p2x4 = np.array([0,1])
p2x7 = np.array([0,1])
p2yc = np.array([0,1])
p2y12 = np.array([0,1])
A = np.array([10,50,100,1000])

time = 300
counter = 0
for i in np.arange(len(A)): # ATP
    for j in np.arange(len(p2yc)): # P2Yc
        for k in np.arange(len(p2x7)): # P2X7
            for l in np.arange(len(p2x4)): # P2X4
                for n in np.arange(len(p2y12)): #P2Y12
                    if p2yc[j] + p2x7[k] + p2x4[l] + p2y12[n] == 0:
                        continue
                    elif p2x7[k] + p2x4[l] + p2y12[n] == 0:
                        continue 
                    elif p2yc[j] + p2x4[l] + p2y12[n] == 0:
                        continue 
                    elif p2yc[j] + p2x7[k] + p2y12[n] == 0:
                        continue
                    elif p2yc[j] + p2x7[k] + p2x4[l]  == 0:
                        continue
                    elif p2x4[l] + p2y12[n]  == 0:
                        continue
                    elif p2yc[j] + p2x7[k]  == 0:
                        continue
                    elif p2yc[j] + p2y12[n]  == 0:
                        continue
                    elif p2yc[j] + p2x4[l] == 0:
                        continue
                    elif p2x7[k] + p2y12[n] == 0:
                        continue
                    elif p2x7[k] + p2x4[l]== 0:
                        continue
                        
                    data   = SR.gotranMicroglia(sim_time      = time,
                                                ATP           = A[i],
                                                output_name   = 'test1',
                                                ode_file_name = 'p2xp2yMigration29', ## it was 22 
                                                data_name2    = 'mRNA_TNF',
                                                data_name3    = 'PI3Ka',
                                                data_name4    = 'Ca4_CN',
                                                data_name5    = 'CaMCN',
                                                data_name6    = 'pAkt',
                                                data_name7    = 'S',
                                                data_name8    = 'NFATpn',
                                                data_name9    = 'NFATNn',
                                                data_name10   = 'TNFae',
                                                rhop2x4       = p2x4[l],
                                                rhop2x7       = p2x7[k],
                                                rhop2yc       = p2yc[j],
                                                rhop2y12      = p2y12[n],
                                                DegSwitch     = 1,
                                                removePickle  = 1)
                    
                    if counter == 0:
                        dura    = data[0]
                        Ca      = data[1]
                        mRNATNF = data[2]
                        PI3K    = data[3]
                        actCN   = 0.1*data[4] + data[5]
                        pakt    = data[6]
                        dist    = data[7]
                        nfatpn  = data[8]
                        nfatnn   = data[9]
                        TNFae    = data[10]
                        entry   = counter
                        atp     = A[i]
                        p4      = p2x4[l]
                        p7      = p2x7[k]
                        py      = p2yc[j]
                        p12     = p2y12[n]
                    else:
                        Ca      = np.vstack([Ca,data[1]])
                        mRNATNF = np.vstack([mRNATNF,data[2]])
                        PI3K    = np.vstack([PI3K,data[3]])
                        actCN   = np.vstack([actCN,0.1*data[4] + data[5]])
                        pakt    = np.vstack([pakt,data[6]])
                        dist    = np.vstack([dist,data[7]])
                        nfatpn  = np.vstack([nfatpn,data[8]])
                        nfatnn   = np.vstack([nfatnn,data[9]])
                        TNFae    = np.vstack([TNFae,data[10]])
                        entry   = np.append(entry,counter)
                        atp     = np.append(atp,A[i])
                        p4      = np.append(p4,p2x4[l])
                        p7      = np.append(p7,p2x7[k])
                        py      = np.append(py,p2yc[j])
                        p12     = np.append(p12,p2y12[n])
                              
                    counter = counter + 1
                    
                    
objects = ('10 uM','50uM','100 uM','1000 uM')
y_pos = np.arange(len(objects))
x = 4 # control
y = 3  # P2Y12 -/-
z = 1  # P2X7 -/-
w = 2  # P2X4 -/-
t = 0  # P2Yc -/-
r = 5
data = np.array([mRNATNF[x][-1],mRNATNF[x+r][-1],mRNATNF[x+r*2][-1],mRNATNF[x+r*3][-1],
                 mRNATNF[y][-1],mRNATNF[y+r][-1],mRNATNF[y+r*2][-1],mRNATNF[y+r*3][-1],
                 mRNATNF[z][-1],mRNATNF[z+r][-1],mRNATNF[z+r*2][-1],mRNATNF[z+r*3][-1],
                 mRNATNF[w][-1],mRNATNF[w+r][-1],mRNATNF[w+r*2][-1],mRNATNF[w+r*3][-1],
                 mRNATNF[t][-1],mRNATNF[t+r][-1],mRNATNF[t+r*2][-1],mRNATNF[t+r*3][-1]])

control = (np.array([mRNATNF[x][-1],mRNATNF[x+r][-1],mRNATNF[x+r*2][-1],mRNATNF[x+r*3][-1]]) - min(data))/(max(data) - min(data))
p2y12 = (np.array([mRNATNF[y][-1],mRNATNF[y+r][-1],mRNATNF[y+r*2][-1],mRNATNF[y+r*3][-1]]) - min(data))/(max(data) - min(data))
p2x7 = (np.array([mRNATNF[z][-1],mRNATNF[z+r][-1],mRNATNF[z+r*2][-1],mRNATNF[z+r*3][-1]]) - min(data))/(max(data) - min(data))
p2x4 = (np.array([mRNATNF[w][-1],mRNATNF[w+r][-1],mRNATNF[w+r*2][-1],mRNATNF[w+r*3][-1]]) - min(data))/(max(data) - min(data))
p2yc = (np.array([mRNATNF[t][-1],mRNATNF[t+r][-1],mRNATNF[t+r*2][-1],mRNATNF[t+r*3][-1]]) - min(data))/(max(data) - min(data))

bar_width = 0.1
space = bar_width/2

plt.figure(figsize=(21,4),dpi=200)
ax = plt.subplot(1,3,1)
plt.title('A',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-space*4, control, bar_width, align='center', color='white', alpha=1, edgecolor='black', label='Control')
plt.bar(y_pos-space*2, p2y12, bar_width, align='center', color='gainsboro', alpha=1, hatch='\\', edgecolor='black', label='P2Y12 -/-')
plt.bar(y_pos+space*0, p2x7, bar_width, align='center', color='gray', alpha=1, edgecolor='black', label='P2X7 -/-')
plt.bar(y_pos+space*2, p2x4, bar_width, align='center', color='black', alpha=1, edgecolor='black', label='P2X4 -/-')
plt.bar(y_pos+space*4, p2yc, bar_width, align='center', color='gainsboro', alpha=1, hatch='//', edgecolor='black', label='P2Yc -/-')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('mRNA TNFa [relative scale]',fontsize=12)
plt.xlabel('[ATP]',fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()

data = np.array([TNFae[x][-1],TNFae[x+r][-1],TNFae[x+r*2][-1],TNFae[x+r*3][-1],
                 TNFae[y][-1],TNFae[y+r][-1],TNFae[y+r*2][-1],TNFae[y+r*3][-1],
                 TNFae[z][-1],TNFae[z+r][-1],TNFae[z+r*2][-1],TNFae[z+r*3][-1],
                 TNFae[w][-1],TNFae[w+r][-1],TNFae[w+r*2][-1],TNFae[w+r*3][-1],
                 TNFae[t][-1],TNFae[t+r][-1],TNFae[t+r*2][-1],TNFae[t+r*3][-1]])
control = (np.array([TNFae[x][-1],TNFae[x+r][-1],TNFae[x+r*2][-1],TNFae[x+r*3][-1]]) - min(data))/(max(data) - min(data))
p2y12 = (np.array([TNFae[y][-1],TNFae[y+r][-1],TNFae[y+r*2][-1],TNFae[y+r*3][-1]]) - min(data))/(max(data) - min(data))
p2x7 = (np.array([TNFae[z][-1],TNFae[z+r][-1],TNFae[z+r*2][-1],TNFae[z+r*3][-1]]) - min(data))/(max(data) - min(data))
p2x4 = (np.array([TNFae[w][-1],TNFae[w+r][-1],TNFae[w+r*2][-1],TNFae[w+r*3][-1]]) - min(data))/(max(data) - min(data))
p2yc = (np.array([TNFae[t][-1],TNFae[t+r][-1],TNFae[t+r*2][-1],TNFae[t+r*3][-1]]) - min(data))/(max(data) - min(data))

ax = plt.subplot(1,3,2)
plt.title('B',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-space*4, control, bar_width, align='center', color='white', alpha=1, edgecolor='black', label='Control')
plt.bar(y_pos-space*2, p2y12, bar_width, align='center', color='gainsboro', alpha=1, hatch='\\', edgecolor='black', label='P2Y12 -/-')
plt.bar(y_pos+space*0, p2x7, bar_width, align='center', color='gray', alpha=1, edgecolor='black', label='P2X7 -/-')
plt.bar(y_pos+space*2, p2x4, bar_width, align='center', color='black', alpha=1, edgecolor='black', label='P2X4 -/-')
plt.bar(y_pos+space*4, p2yc, bar_width, align='center', color='gainsboro', alpha=1, hatch='//', edgecolor='black', label='P2Yc -/-')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('Released TNFa [relative scale]',fontsize=12)
plt.xlabel('[ATP]',fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()

data = np.array([dist[x][-1],dist[x+r][-1],dist[x+r*2][-1],dist[x+r*3][-1],
                 dist[y][-1],dist[y+r][-1],dist[y+r*2][-1],dist[y+r*3][-1],
                 dist[z][-1],dist[z+r][-1],dist[z+r*2][-1],dist[z+r*3][-1],
                 dist[w][-1],dist[w+r][-1],dist[w+r*2][-1],dist[w+r*3][-1],
                 dist[t][-1],dist[t+r][-1],dist[t+r*2][-1],dist[t+r*3][-1]])

#control = (np.array([dist[x][-1],dist[x+r][-1],dist[x+r*2][-1],dist[x+r*3][-1]]) - min(data))/(max(data) - min(data))
#p2y12 = (np.array([dist[y][-1],dist[y+r][-1],dist[y+r*2][-1],dist[y+r*3][-1]]) - min(data))/(max(data) - min(data))
#p2x7 = (np.array([dist[z][-1],dist[z+r][-1],dist[z+r*2][-1],dist[z+r*3][-1]]) - min(data))/(max(data) - min(data))
#p2x4 = (np.array([dist[w][-1],dist[w+r][-1],dist[w+r*2][-1],dist[w+r*3][-1]]) - min(data))/(max(data) - min(data))
#p2xc = (np.array([dist[t][-1],dist[t+r][-1],dist[t+r*2][-1],dist[t+r*3][-1]]) - min(data))/(max(data) - min(data))
control = np.array([dist[x][-1],dist[x+r][-1],dist[x+r*2][-1],dist[x+r*3][-1]])/(max(data))
p2y12   = np.array([dist[y][-1],dist[y+r][-1],dist[y+r*2][-1],dist[y+r*3][-1]])/(max(data))
p2x7    = np.array([dist[z][-1],dist[z+r][-1],dist[z+r*2][-1],dist[z+r*3][-1]])/(max(data))
p2x4    = np.array([dist[w][-1],dist[w+r][-1],dist[w+r*2][-1],dist[w+r*3][-1]])/(max(data))
p2yc    = np.array([dist[t][-1],dist[t+r][-1],dist[t+r*2][-1],dist[t+r*3][-1]])/(max(data))


ax = plt.subplot(1,3,3)
plt.title('C',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-space*4, control, bar_width, align='center', color='white', alpha=1, edgecolor='black', label='Control')
plt.bar(y_pos-space*2, p2y12, bar_width, align='center', color='gainsboro', alpha=1, hatch='\\', edgecolor='black', label='P2Y12 -/-')
plt.bar(y_pos+space*0, p2x7, bar_width, align='center', color='gray', alpha=1, edgecolor='black', label='P2X7 -/-')
plt.bar(y_pos+space*2, p2x4, bar_width, align='center', color='black', alpha=1, edgecolor='black', label='P2X4 -/-')
plt.bar(y_pos+space*4, p2yc, bar_width, align='center', color='gainsboro', alpha=1, hatch='//', edgecolor='black', label='P2Yc -/-')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('Scale of Chemotaxis [a.u.]',fontsize=12)
plt.xlabel('ATP',fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()

plt.subplots_adjust(wspace=0.25)
plt.savefig("mRNAvsTNFvsMig.png")

## SA test for elevation vs frequency

time = 600
counter = 0
Kp = scipy.linspace(5,20,10)
p2x7 = scipy.linspace(0.01,1,10)
for i in np.arange(len(Kp)): # ATP
    data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = 200,
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration29', # 6 works
                                    DegSwitch     = 0,
                                    kg_p2y        = Kp[i],
                                    data_name2    = 'TNFae',
                                    data_name3    = 'mRNA_TNF',
                                    data_name4    = 'S',
                                    removePickle  = 1)
    
    if counter == 0:
            dura    = data[0]
            Ca      = data[1]
            Catail  = np.average(data[1][300:599])
            tnfe    = data[2]
            mrna    = data[3]
            dist    = data[4]
            entry   = counter
            kp      = Kp[i]
    else:
            Ca      = np.vstack([Ca,data[1]])
            Catail  = np.append(Catail,np.average(data[1][300:599]))
            tnfe    = np.vstack([tnfe,data[2]])
            mrna    = np.vstack([mrna,data[3]])
            dist    = np.vstack([dist,data[4]])
            kp      = np.append(kp,Kp[i])
        
    counter = counter + 1 

distdummy = np.zeros(10)
tnfdummy = np.zeros(10)
for i in np.arange(10):
    distdummy[i] = dist[i][300]
    tnfdummy[i] = tnfe[i][300]
    

time = 600
counter = 0
p2x7 = scipy.linspace(0.01,50,10)
for j in np.arange(len(p2x7)):
    data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = 200,
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration29', # 6 works
                                    rhop2x7       = p2x7[j],
                                    DegSwitch     = 0,
                                    EleSwitch     = 1,
                                    rhop2yc       = 0,
                                    data_name2    = 'TNFae',
                                    data_name3    = 'mRNA_TNF',
                                    data_name4    = 'S',
                                    data_name5    = 'NFATNn',
                                    removePickle  = 1)
    
    if counter == 0:
            dura    = data[0]
            Ca      = data[1]
            Catail  = np.average(data[1][300:599])
            tnfe    = data[2]
            mrna    = data[3]
            dist    = data[4]
            nfat    = data[5]
            entry   = counter
            p7      = p2x7[j]
    else:
            Ca      = np.vstack([Ca,data[1]])
            Catail  = np.append(Catail,np.average(data[1][300:599]))
            tnfe    = np.vstack([tnfe,data[2]])
            mrna    = np.vstack([mrna,data[3]])
            dist    = np.vstack([dist,data[4]])
            nfat    = np.vstack([nfat,data[5]])
            p7      = np.append(p7,p2x7[j])
        
    counter = counter + 1 
    
    
distdummy2 = np.zeros(10)
tnfdummy2 = np.zeros(10)
nfatdummy2 = np.zeros(10)
for i in np.arange(10):
    distdummy2[i] = dist[i][300]
    tnfdummy2[i] = tnfe[i][300]
    nfatdummy2[i] = nfat[i][300]
    
datadummy = distdummy-min(distdummy)
datas1 = datadummy/max(datadummy)
datadummy = distdummy2-min(distdummy2)
datas2 = datadummy/max(datadummy)
freq = scipy.linspace(0,2.7,10)
plt.figure(figsize=(6,4),dpi=100)
plt.tick_params(direction='in',labelsize=12)
plt.plot(freq,datas1,'bs-',lw=2,alpha=1)
plt.xlabel('Freq. [peak per min]',fontsize=15)
plt.ylabel('Scale of Chemotaxis [a.u.]',fontsize=15)
plt.tight_layout()
plt.savefig("chemfreq.png")

plt.figure(figsize=(6,4),dpi=100)
plt.tick_params(direction='in',labelsize=12)
plt.plot(Catail,datas2,'bs-',lw=2,alpha=1)
plt.xlabel('$Ca^{2+}$ Elevation [nM]',fontsize=15)
plt.ylabel('Scale of Chemotaxis [a.u.]',fontsize=15)
plt.tight_layout()
plt.savefig("chemelev.png")

datadummy = tnfdummy-min(tnfdummy)
datas3 = datadummy/max(datadummy)
datadummy = tnfdummy2-min(tnfdummy2)
datas4 = datadummy/max(datadummy)
datadummy = nfatdummy2-min(nfatdummy2)
datas5 = datadummy/datadummy[-2]
freq = scipy.linspace(0,2.7,10)

plt.figure(figsize=(6,4),dpi=100)
plt.tick_params(direction='in',labelsize=12)
plt.plot(freq,datas3,'bs-',lw=2,alpha=1)
plt.xlabel('Freq. [peak per min]',fontsize=15)
plt.ylabel('Released TNFa [scaled by max]',fontsize=15)
plt.tight_layout()
plt.savefig("tnffreq.png")

plt.figure(figsize=(6,4),dpi=100)
plt.tick_params(direction='in',labelsize=12)
plt.plot(Catail,datas4,'bs-',lw=2,alpha=1)
plt.xlabel('$Ca^{2+}$ Elevation [nM]',fontsize=15)
plt.ylabel('Released TNFa [scaled by max]',fontsize=15)
plt.tight_layout()
plt.savefig("tnfelev.png")


### Color Map

time = 600
counter = 0
Kp = scipy.linspace(5,20,20)
p2x7 = scipy.linspace(0.01,1,20)
for i in np.arange(len(Kp)): # ATP
    for j in np.arange(len(p2x7)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = 200,
                                    output_name   = 'test1',
                                    ode_file_name = 'p2xp2yMigration29', # 6 works
                                    rhop2x7       = p2x7[j],
                                    DegSwitch     = 0,
                                    EleSwitch     = 1,
                                    kg_p2y        = Kp[i],
                                    data_name2    = 'TNFae',
                                    data_name3    = 'mRNA_TNF',
                                    data_name4    = 'S',
                                    removePickle  = 1)
    
        if counter == 0:
            dura    = data[0]
            Ca      = data[1]
            Catail  = np.average(data[1][300:599])
            tnfe    = data[2]
            mrna    = data[3]
            dist    = data[4]
            entry   = counter
            p7      = p2x7[j]
            kp      = Kp[i]
        else:
            Ca      = np.vstack([Ca,data[1]])
            Catail  = np.append(Catail,np.average(data[1][300:599]))
            tnfe    = np.vstack([tnfe,data[2]])
            mrna    = np.vstack([mrna,data[3]])
            dist    = np.vstack([dist,data[4]])
            p7      = np.append(p7,p2x7[j])
            kp      = np.append(kp,Kp[i])
        
        counter = counter + 1 
        
        
        
minCa = min(Catail)
maxCa = max(Catail)
Calim = scipy.linspace(minCa,maxCa,30)
m = len(Calim)
w = len(Kp)
l = len(Catail)
distA = np.zeros((m,w))
tnfeA = np.zeros((m,w))
mrnaA = np.zeros((m,w))
Calist = np.zeros((m,w))
Kplist = np.zeros((m,w))
counter = 0
for i in np.arange(w): ## frequency base
    for j in np.arange(m): ## concentration base
        for k in  np.arange(l): ## data length base
            if Kp[(w-1)-i] == kp[k]:
                if j+1 == m:
                    continue
                else:
                    if Catail[k] > Calim[j] and Catail[k] < Calim[j+1]:
                        
                        distA[j,i] = dist[k][299]
                        tnfeA[j,i] = tnfe[k][299]
                        mrnaA[j,i] = mrna[k][299]
                        Kplist[j,i] = Kp[w-1-i]
                        Calist[j,i] = Catail[k]
                    
for i in np.arange(w):
    for j in np.arange(m):
        if distA[j,i] == 0:
            distA[j,i] = float('Nan')
            tnfeA[j,i] = float('Nan')
            mrnaA[j,i] = float('Nan')
            counter = counter + 1
            
            
import matplotlib.cm as cm
plt.figure(figsize=(7,5),dpi=150)
plt.tick_params(labelsize=12,direction='in')
## Control 3 = 1/5 1 = 1/30 2 = 1/60
nan = 0 #float('NaN')

data2 = (np.array(np.flip(distA,1))-distA[0][-1])/(distA[28][1]-distA[0][-1])

data = scipy.ndimage.zoom(data2, 4)

im = plt.imshow(data2,  origin='lower',
                cmap=cm.winter, extent=(0,2.7,minCa,maxCa),aspect='auto')
#im = plt.imshow(data2, origin='lower',
#                cmap=cm.winter, extent=(1,4,minCa,maxCa),aspect='auto')

#levels = [0.5]
#CS = plt.contour(data2, levels, origin='lower',colors='yellow', 
#                linewidths=2, extent=(1,4,minCa,maxCa))

CBI = plt.colorbar(im, orientation='vertical', shrink=0.8)
CBI.set_label('Degree of Chemotaxis [scale]', fontsize=12)
#plt.clabel(CS, inline=1, fontsize=12)
#plt.xlim([0,3])
#plt.ylim([100,106])
plt.xlabel("Frequency [peak per min]",fontsize=15)
plt.ylabel("$[Ca^{2+}]_i$ elevation [nM]",fontsize=15)
plt.tight_layout()
plt.savefig('freqvselevDist.png')

import matplotlib.cm as cm
plt.figure(figsize=(7,5),dpi=150)
plt.tick_params(labelsize=12,direction='in')
## Control 3 = 1/5 1 = 1/30 2 = 1/60
nan = 0 #float('NaN')

data2 = (np.array(np.flip(tnfeA,1))-tnfeA[0][-1])/(tnfeA[28][1]-tnfeA[0][-1])

#data = scipy.ndimage.zoom(data2, 4)

im = plt.imshow(data2,  origin='lower',
                cmap=cm.winter, extent=(0,2.7,minCa,Calim[9]),aspect='auto')
#im = plt.imshow(data2, origin='lower',
#                cmap=cm.winter, extent=(1,4,minCa,maxCa),aspect='auto')

#levels = [0.5]
#CS = plt.contour(data2, levels, origin='lower',colors='yellow', 
#                linewidths=2, extent=(1,4,minCa,maxCa))

CBI = plt.colorbar(im, orientation='vertical', shrink=0.8)
CBI.set_label('Normalized Released TNFa [scale]', fontsize=12)
#plt.clabel(CS, inline=1, fontsize=12)
#plt.xlim([0,3])
plt.xlabel("Frequency [peak per min]",fontsize=15)
plt.ylabel("$[Ca^{2+}]_i$ elevation [nM]",fontsize=15)
plt.tight_layout()
plt.savefig('freqvselevTNFaReleased.png')

# Primary Cultured
dataPrime10   = SR.gotranMicroglia(sim_time      = 300,
                                 ATP           = 10,
                                 rhop2yc       = 12/15,
                                 rhop2y12      = 0,
                                 rhop2x4       = 3,
                                 rhop2x7       = 0.55,
                                 output_name   = 'test1',
                                 ode_file_name = 'p2xp2yMigration29',
                                 data_name2    = 'pAkt',
                                 data_name3    = 'S',
                                 data_name4    = 'TNFae',
                                 removePickle  = 1,
                                 timePrint     = 0)
dataPrime100   = SR.gotranMicroglia(sim_time      = 300,
                                 ATP           = 100,
                                 rhop2yc       = 12/15,
                                 rhop2y12      = 0,
                                 rhop2x4       = 3,
                                 rhop2x7       = 0.55,
                                 output_name   = 'test1',
                                 ode_file_name = 'p2xp2yMigration29',
                                 data_name2    = 'pAkt',
                                 data_name3    = 'S',
                                 data_name4    = 'TNFae',
                                 removePickle  = 1,
                                 timePrint     = 0)
dataPrime1000   = SR.gotranMicroglia(sim_time      = 300,
                                 ATP           = 1000,
                                 rhop2yc       = 12/15,
                                 rhop2y12      = 0,
                                 rhop2x4       = 3,
                                 rhop2x7       = 0.55,
                                 output_name   = 'test1',
                                 ode_file_name = 'p2xp2yMigration29',
                                 data_name2    = 'pAkt',
                                 data_name3    = 'S',
                                 data_name4    = 'TNFae',
                                 removePickle  = 1,
                                 timePrint     = 0)

# Resting
dataRest10   = SR.gotranMicroglia(sim_time      = 300,
                                ATP           = 10,
                                rhop2yc       = 1,
                                rhop2y12      = 1,
                                rhop2x4       = 1,
                                rhop2x7       = 1,
                                output_name   = 'test1',
                                ode_file_name = 'p2xp2yMigration29',
                                data_name2    = 'pAkt',
                                data_name3    = 'S',
                                data_name4    = 'TNFae',
                                removePickle  = 1,
                                timePrint     = 0)
dataRest100   = SR.gotranMicroglia(sim_time      = 300,
                                ATP           = 100,
                                rhop2yc       = 1,
                                rhop2y12      = 1,
                                rhop2x4       = 1,
                                rhop2x7       = 1,
                                output_name   = 'test1',
                                ode_file_name = 'p2xp2yMigration29',
                                data_name2    = 'pAkt',
                                data_name3    = 'S',
                                data_name4    = 'TNFae',
                                removePickle  = 1,
                                timePrint     = 0)
dataRest1000   = SR.gotranMicroglia(sim_time      = 300,
                                ATP           = 1000,
                                rhop2yc       = 1,
                                rhop2y12      = 1,
                                rhop2x4       = 1,
                                rhop2x7       = 1,
                                output_name   = 'test1',
                                ode_file_name = 'p2xp2yMigration29',
                                data_name2    = 'pAkt',
                                data_name3    = 'S',
                                data_name4    = 'TNFae',
                                removePickle  = 1,
                                timePrint     = 0)


pakt10 = np.array([dataRest10[2][-1],dataPrime10[2][-1]])/dataRest1000[2][-1]
pakt100 = np.array([dataRest100[2][-1],dataPrime100[2][-1]])/dataRest1000[2][-1]
pakt1000 = np.array([dataRest1000[2][-1],dataPrime1000[2][-1]])/dataRest1000[2][-1]
migr10 = np.array([dataRest10[3][-1],dataPrime10[3][-1]])/dataRest1000[3][-1]
migr100 = np.array([dataRest100[3][-1],dataPrime100[3][-1]])/dataRest1000[3][-1]
migr1000 = np.array([dataRest1000[3][-1],dataPrime1000[3][-1]])/dataRest1000[3][-1]
tnfa10 = (np.array([dataRest10[4][-1],dataPrime10[4][-1]])-min(dataRest10[4]))/(dataRest1000[4][-1]-min(dataRest1000[4]))
tnfa100 = (np.array([dataRest100[4][-1],dataPrime100[4][-1]])-min(dataRest10[4]))/(dataRest1000[4][-1]-min(dataRest1000[4]))
tnfa1000 = (np.array([dataRest1000[4][-1],dataPrime1000[4][-1]])-min(dataRest10[4]))/(dataRest1000[4][-1]-min(dataRest1000[4]))

objects = ('Resting','Primary')
y_pos = np.arange(len(objects))

bar_width = 0.25
plt.figure(figsize=(8,4),dpi=150)
ax = plt.subplot(1,2,1)
plt.title('A',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-0.25, tnfa10, bar_width, align='center', color='black', alpha=1, edgecolor='black',label='10 uM')
plt.bar(y_pos,      tnfa100, bar_width, align='center', color='gray', alpha=1, edgecolor='black',label='100 uM')
plt.bar(y_pos+0.25, tnfa1000, bar_width, align='center', color='white', alpha=1, edgecolor='black',label='1000 uM')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('Released TNFa [scale to max]',fontsize=12)
plt.tight_layout()

ax = plt.subplot(1,2,2)
plt.title('B',fontsize=12,fontweight='bold',loc='left')
plt.bar(y_pos-0.25, migr10, bar_width, align='center', color='black', alpha=1, edgecolor='black',label='10 uM')
plt.bar(y_pos,      migr100, bar_width, align='center', color='gray', alpha=1, edgecolor='black',label='100 uM')
plt.bar(y_pos+0.25, migr1000, bar_width, align='center', color='white', alpha=1, edgecolor='black',label='1000 uM')
plt.tick_params(labelsize=12,direction='in')
plt.xticks(y_pos, objects)
plt.ylabel('Mig. Dist./Chemotaxis [scale to max]',fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.tight_layout()

plt.subplots_adjust(wspace=0.25)
plt.savefig("RestvsPrime.png")

