#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 11:12:42 2025

@author: arthurarnoult
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import os
os.environ['PATH'] += '/Library/TeX/texbin'

rc('text', usetex=True)
rc('font', family='serif')

I=5000
K=50



def Histogramme(li,li_nom_label):
#    data=np.zeros((I,len(li_nom_fichier)))
#    
#    for i in range(len(li_nom_fichier)):
#        data[:,i]=Open_Numbers_file(li_nom_fichier[i])
    
    data=li
 
    
    counts, bins, _ = plt.hist(data, bins=100, density=True, histtype='step', label=li_nom_label)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1])
    
    plt.xlabel("Répartition des $\mathcal{G}$")
    plt.ylabel("Densité")


    plt.title("$T=50$, $N=1000$, $X_0=1$, $I=25000$")
    
    plt.savefig("Densité_poids_log.pdf")
    
    
li=np.zeros((I,K+1))
fichier=open("/Users/arthurarnoult/EMC2/Tony/Code_IS_Parareal/RES_1D/DENSITE/G_cal/MC_Parareal_G_T=50_N=1000_K=100_I=25000_Nb_atoms=1_1.txt").readlines()
for i in range(I-1):
   a=fichier[i].split("\t")
   for j in range(K+1):
#       li[i,j]=np.log(float(a[j]))
       li[i,j]=float(a[j])

              

li_densite=np.linspace(800,1200,301)
def f(x):
    return 1/np.sqrt(4000*np.pi)*np.exp(-(x-1000)**2/4000)

plt.plot(li_densite, f(li_densite), label=r'$\mathcal{N}(1000,2000)$')

plt.gca().set_xlim([800,1300])



li_0=li[:,0]
li_5=li[:,5]
li_10=li[:,10]
li_20=li[:,20]
li_50=li[:,50]
#li_100=li[:,100]

li_extrait=[li_0,li_5,li_10,li_20, li_50]
Histogramme(li_extrait,['$k=0$','$k=5$','$k=10$','$k=20$','$k=50$'])
plt.savefig("Analyse/FIG/Densité_Gcal.pdf")


#Histogramme(li_5,['$k=5$'])