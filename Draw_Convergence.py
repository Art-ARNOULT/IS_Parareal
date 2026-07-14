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

I=10000
K=20

#fichier_esp=open("RES/Esp_T=5_N=100_X0=1.5_K=50_I="+str(I)+".txt").readlines()
#fichier_IC=open("RES/IC_T=5_N=100_X0=1.5_K=50_I="+str(I)+".txt").readlines()

fichier_esp=open("RES/Esp_T=5_N=100_X0=1.5_K=20_I=10000_weight_sigma=x.txt").readlines()
fichier_re=open("RES/IC_T=5_N=100_X0=1.5_K=20_I=10000_weight_sigma=x.txt").readlines()

fichier_IC=open("RES/IC_T=5_N=100_X0=1.5_K=20_I=10000_weight_sigma=x.txt").readlines()


I=len(fichier_esp)


esp=np.zeros((I,K+1))
esp_re=np.zeros((I,K+1))
IC=np.zeros((I,K+1))


Valeur_I=np.zeros(I)
for i in range(I):
   a=fichier_esp[i].split("\t")
   b=fichier_IC[i].split("\t")
#   c=fichier_re[i].split("\t")
   Valeur_I[i]=a[0]
   
   for k in range(1,K+2):
       esp[i,k-1]=float(a[k])
       IC[i,k-1] =float(b[k])
#       esp_re[i,k-1]=float(c[k])



#fichier_esp_ecriture=open("Esp_T=130_N=2600_X0=1_K=100_I=100000_PSI_save.txt", "w+")
#fichier_esp_re_ecriture=open("Esp_reweighting_T=130_N=2600_X0=1_K=100_I=100000_PSI_save.txt", "w+")
#
#for k in range(K):
#    fichier_esp_ecriture.write(str(k) + "\t" + str(esp[-1,k])+"\n")
#    fichier_esp_re_ecriture.write(str(k) + "\t" + str(esp_re[-1,k])+"\n")
#fichier_esp_ecriture.close()
#fichier_esp_re_ecriture.close()



#x=np.arange(0,K+1)
##plt.plot(esp[-1,:], color='green', label="k="+str(k))
#plt.title("$T=5$, $N=100$, $X_0=1.5$")
#plt.plot(x,esp[-1,:], color='red', label="Normal")
##plt.plot(x,esp_re[-1,:], color='green', label="reweighting")
#plt.legend()
#plt.yscale("log")
#plt.axhline(4.22, color="black", alpha=0.2,linestyle="--")
#plt.savefig("FIG/CV_selon_k_psi_T=130_N=2300_X0=1_I="+str(I)+".pdf", dpi=200)
#plt.gca().set_ylim([0,50])
#
#plt.show()



for k in range(0,K+1,5):
    
    fig=plt.figure()
    plt.title("$T=5$, $N=100$, $X_0=1.5$")
    plt.xlabel("Itération Monte-Carlo I")
    plt.ylabel("Espérence empirique")
    

#    x_decim = x[::100]
#    esp_decim = esp[:,k][::100]
    # IC_decim=IC[:,k][::100]

#    plt.gca().set_ylim([0.5,1.5])

    plt.plot(Valeur_I, esp[:,k], color='green', label="k="+str(k))
    plt.fill_between(Valeur_I,  esp[:,k]- IC[:,k], esp[:,k] + IC[:,k], color='green', alpha=0.2)
    plt.title("$T=5$, $N=100$, $X_0=1.5$ -- Poids, $\sigma(x)= x$")
    plt.legend()
    plt.axhline(1, color="black", alpha=0.2,linestyle="--")
    # plt.savefig("FIG/CV_MC_poids_T=5_N=100_X0=1.5_k="+str(k)+"_I="+str(I)+".png", dpi=200)
#    plt.show()
    plt.savefig("CV_poids_sigma=x_k="+str(k)+".pdf")
