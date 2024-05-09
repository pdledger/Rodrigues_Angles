# Python script for computing rodrigues angles and comparing different angle measures

# Routine: main.py
# Author: P.D. Ledger
# Date: 18/1/2024

# include a few standard libaries
from time import time
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

# include the src path
sys.path.insert(0, "src")


from EigenValueDecomposition import *
from MinMaxthetafromQRQI import *
from SortEigenValues import *
from AnglesSortedQRQI import *
from SplitTensor import *
from Fmeasure import *
from Commeasure import *



#   Main script pass arguments of
#   directory = Directory in which Tensors.csv, N0.csv and Frequencies.csv are stored
#   MaxOmega = Maximum frequency to be considered
def main_noneigveconly(directory,MaxOmega,Figures="On",FullRom="Rom"):
    print("Opening files from this path = ",directory)

    # Read data
    if FullRom=="Full":
        TensorArray = np.genfromtxt(directory+'PODTensors.csv', dtype=complex, delimiter=', ')
        N0 =          np.genfromtxt(directory+'N0.csv', dtype=float, delimiter=',')
        Frequencies = np.genfromtxt(directory+'PODFrequencies.csv', dtype=float, delimiter=', ')
    elif FullRom=="Rom":
        TensorArray = np.genfromtxt(directory+'Tensors.csv', dtype=complex, delimiter=', ')
        N0 =          np.genfromtxt(directory+'N0.csv', dtype=float, delimiter=',')
        Frequencies = np.genfromtxt(directory+'Frequencies.csv', dtype=float, delimiter=', ')

    # Limit the values up to a perscribed MaxOmega
    N=len(Frequencies)
    Ntarget=N
    for n in range(N):
        if Frequencies[n] > MaxOmega:
            Ntarget=n
            break

    # allow for frequency array less than max omega.
    if MaxOmega > np.max(Frequencies):
        Ntarget = N

    TensorArray = TensorArray[:Ntarget,:]
    Frequencies = Frequencies[:Ntarget]
    N=Ntarget

    Rstore,Istore,Rtildestore, N0store = SplitTensor(TensorArray,Frequencies,N0)

    # Determine eigenvalue decompositions of N0, R, I, Rtilde (no sorting applied), and their multiplicities
    MultRstore, MultIstore, MultRtildestore, MultN0store, URstore, UIstore, URtildestore, UN0store, QRstore, QIstore, QRtildestore, QN0store = EigenValueDecomposition(N0,TensorArray,Frequencies)

    # Determine the maximal and minimal angles from QR and QI also output d_F metric for these orderings
    FixEvecs="Yes"#"Yes"
    #MinAnglestoreRI, MaxAnglestoreRI, dFMinAnglestoreRI, dFMaxAnglestoreRI,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = MinMaxthetafromQRQI(Frequencies,QRstore,QIstore,URstore, UIstore,MultRstore,MultIstore,FixEvecs)
    SortedMultRstore=[]
    SortedMultIstore=[]

    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues, Rstore, Istore)

    # Obtain angles for this sorted min-max combination
    #AnglestoreRIsortedmaxdiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreRIfmeasfullconstsortedmaxdiff, AnglestoreRIfmeasapprxconstsortedmaxdiff_min,AnglestoreRIfmeasapprxconstsortedmaxdiff_max, RIfmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    #Obtain Com-meauses (approx and exact constant)
    AnglestoreRIcommeasfullconstsortedmaxdiff, AnglestoreRIcommeasapprxconstsortedmaxdiff_min,AnglestoreRIcommeasapprxconstsortedmaxdiff_max, RIcommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

        # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues, Rstore, Istore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreRIsortedmindiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreRIfmeasfullconstsortedmindiff, AnglestoreRIfmeasapprxconstsortedmindiff_min, AnglestoreRIfmeasapprxconstsortedmindiff_max, RIfmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    #Obtain Com-meauses (approx and exact constant)
    AnglestoreRIcommeasfullconstsortedmindiff, AnglestoreRIcommeasapprxconstsortedmindiff_min,AnglestoreRIcommeasapprxconstsortedmindiff_max, RIcommeapprx_den_const_min = Commeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    kthetavec=np.zeros((N,3))
    #for n in range(N):
    #    kthetavec[n,0] = -SortedKstore[n,1,2]*MinAnglestoreRI[n]
    #    kthetavec[n,1] = SortedKstore[n,0,2]*MinAnglestoreRI[n]
    #    kthetavec[n,2] = -SortedKstore[n,0,1]*MinAnglestoreRI[n]
    if Figures=="On":
        fig=plt.figure()
        #plt.semilogx(Frequencies,MinAnglestoreRI,label=r'$d_R({\cal R},{\cal I})$')
        #plt.semilogx(Frequencies,dFMinAnglestoreRI,label=r'$d_F({\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRIfmeasapprxconstsortedmaxdiff_min,AnglestoreRIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_E({\cal R},{\cal I})$ ')

        #plt.semilogx(Frequencies,np.fmin(AnglestoreRIfmeasapprxconstsortedmaxdiff_max,AnglestoreRIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_E({\cal R},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRIcommeasapprxconstsortedmaxdiff_min,AnglestoreRIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_C({\cal R},{\cal I})$ ')

        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()


        # Plot of vairation of |theta| kvec with omgea
        #fig=plt.figure()
        #plt.semilogx(Frequencies,abs(kthetavec[:,0]),label=r'$|\theta k_1|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,1]),label=r'$|\theta k_2|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,2]),label=r'$|\theta k_3|$')

        #plt.xlabel(r'$\omega$ [rad/s]')
        #plt.ylabel(r'$|\theta k_i| $ [rad]')
        #plt.legend()
        #plt.show()


    RIResults={"Frequencies": Frequencies,
        "AnglestoreRIcommeasapprxconstsortedmaxdiff_max":AnglestoreRIcommeasapprxconstsortedmaxdiff_max,"AnglestoreRIcommeasapprxconstsortedmindiff_min":AnglestoreRIcommeasapprxconstsortedmindiff_min,\
        "AnglestoreRIcommeasapprxconstsortedmaxdiff_min":AnglestoreRIcommeasapprxconstsortedmaxdiff_min,"AnglestoreRIcommeasapprxconstsortedmindiff_max":AnglestoreRIcommeasapprxconstsortedmindiff_max, \
        "RIfmeasapprx_den_const_max": RIfmeasapprx_den_const_max, "RIfmeasapprx_den_const_min": RIfmeasapprx_den_const_min, \
        "RIcommeapprx_den_const_max":RIcommeapprx_den_const_max, "RIcommeapprx_den_const_min": RIcommeapprx_den_const_min, \
        "URstore":URstore, "UIstore":UIstore, "URtildestore":URtildestore, "UN0store":UN0store, "QRstore":QRstore, "QIstore":QIstore, "QRtildestore":QRtildestore, \
         "QN0store":QN0store,"AnglestoreRIfmeasfullconstsortedmindiff": AnglestoreRIfmeasfullconstsortedmindiff,\
         "AnglestoreRIfmeasfullconstsortedmaxdiff": AnglestoreRIfmeasfullconstsortedmaxdiff,"AnglestoreRIcommeasfullconstsortedmindiff":AnglestoreRIcommeasfullconstsortedmindiff,\
         "AnglestoreRIcommeasfullconstsortedmaxdiff":AnglestoreRIcommeasfullconstsortedmaxdiff,"Rstore":Rstore,"Istore":Istore,"Rtildestore":Rtildestore, "N0store":N0store,\
         "RIkthetavec":kthetavec}



    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    #MinAnglestoreRtildeI, MaxAnglestoreRtildeI, dFMinAnglestoreRtildeI, dFMaxAnglestoreRtildeI,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore= MinMaxthetafromQRQI(Frequencies,QRtildestore,QIstore,URtildestore, UIstore,MultRtildestore,MultIstore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues, Rtildestore, Istore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreRtildeIsortedmaxdiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    print('Computing F measure Tilde')
    AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max, RtildeIfmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreRtildeIcommeasfullconstsortedmaxdiff, AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min, AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max, RtildeIcommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues, Rtildestore, Istore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreRtildeIsortedmindiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max, RtildeIfmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    AnglestoreRtildeIcommeasfullconstsortedmindiff, AnglestoreRtildeIcommeasapprxconstsortedmindiff_min,AnglestoreRtildeIcommeasapprxconstsortedmindiff_max, RtildeIcommeapprx_den_const_min = Commeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    kthetavec=np.zeros((N,3))
    #for n in range(N):
    #    kthetavec[n,0] = -SortedKstore[n,1,2]*MinAnglestoreRtildeI[n]
    #    kthetavec[n,1] = SortedKstore[n,0,2]*MinAnglestoreRtildeI[n]
    #    kthetavec[n,2] = -SortedKstore[n,0,1]*MinAnglestoreRtildeI[n]

    if Figures=="On":
        fig=plt.figure()
        #plt.semilogx(Frequencies,MinAnglestoreRtildeI,label=r'$d_R(\tilde{\cal R},{\cal I})$')
        #plt.semilogx(Frequencies,dFMinAnglestoreRtildeI,label=r'$d_F(\tilde{\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_E(\tilde{\cal R},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min,AnglestoreRtildeIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_C(\tilde{\cal R},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

        # Plot of vairation of |theta| kvec with omgea
        #fig=plt.figure()
        #plt.semilogx(Frequencies,abs(kthetavec[:,0]),label=r'$|\theta k_1 |$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,1]),label=r'$|\theta k_2|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,2]),label=r'$|\theta k_3|$')

        #plt.xlabel(r'$\omega$ [rad/s]')
        #plt.ylabel(r'$|\theta k_i |$ [rad]')
        #plt.legend()
        #plt.show()

    RtildeIResults={"Frequencies": Frequencies, "AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max": AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreRtildeIfmeasapprxconstsortedmindiff_min": AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,"AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min": AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,
        "AnglestoreRtildeIfmeasapprxconstsortedmindiff_max": AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,
        "AnglestoreRtildeIfmeasfullconstsortedmindiff": AnglestoreRtildeIfmeasfullconstsortedmindiff,\
        "AnglestoreRtildeIfmeasfullconstsortedmaxdiff": AnglestoreRtildeIfmeasfullconstsortedmaxdiff,\
        "AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max":AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max,"AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min":AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min,\
        "AnglestoreRtildeIcommeasapprxconstsortedmindiff_min":AnglestoreRtildeIcommeasapprxconstsortedmindiff_min, "AnglestoreRtildeIcommeasapprxconstsortedmindiff_max":AnglestoreRtildeIcommeasapprxconstsortedmindiff_max, \
        "RtildeIfmeasapprx_den_const_max": RtildeIfmeasapprx_den_const_max, "RtildeIfmeasapprx_den_const_min": RtildeIfmeasapprx_den_const_min,
        "RtildeIcommeapprx_den_const_max":RtildeIcommeapprx_den_const_max, "RtildeIcommeapprx_den_const_min": RtildeIcommeapprx_den_const_min,\
        "RtildeIkthetavec":kthetavec}


    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    #MinAnglestoreN0I, MaxAnglestoreN0I, dFMinAnglestoreN0I, dFMaxAnglestoreN0I,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore = MinMaxthetafromQRQI(Frequencies,QN0store,QIstore,UN0store, UIstore,MultN0store,MultIstore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultN0store, SortedMultIstore, SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore = SortEigenValues(MultN0store, MultIstore, UN0store, UIstore, QN0store, QIstore, Frequencies, sorteigenvalues, N0store, Istore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreN0Isortedmaxdiff = AnglesSortedQRQI(SortedQN0store,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Ifmeasfullconstsortedmaxdiff, AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min, AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max, N0Ifmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreN0Icommeasfullconstsortedmaxdiff, AnglestoreN0Icommeasapprxconstsortedmaxdiff_min, AnglestoreN0Icommeasapprxconstsortedmaxdiff_max, N0Icommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    SortedMultN0store, SortedMultIstore, SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore = SortEigenValues(MultN0store, MultIstore, UN0store, UIstore, QN0store, QIstore, Frequencies, sorteigenvalues, N0store, Istore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreN0Isortedmindiff = AnglesSortedQRQI(SortedQN0store,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Ifmeasfullconstsortedmindiff, AnglestoreN0Ifmeasapprxconstsortedmindiff_min,AnglestoreN0Ifmeasapprxconstsortedmindiff_max, N0Ifmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    AnglestoreN0Icommeasfullconstsortedmindiff, AnglestoreN0Icommeasapprxconstsortedmindiff_min,AnglestoreN0Icommeasapprxconstsortedmindiff_max, N0Icommeapprx_den_const_min= Commeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)

    #kthetavec=np.zeros((N,3))
    #for n in range(N):
    #    kthetavec[n,0] = -SortedKstore[n,1,2]*MinAnglestoreN0I[n]
    #    kthetavec[n,1] = SortedKstore[n,0,2]*MinAnglestoreN0I[n]
    #    kthetavec[n,2] = -SortedKstore[n,0,1]*MinAnglestoreN0I[n]

    if Figures=="On":
        fig=plt.figure()
        #plt.semilogx(Frequencies,MinAnglestoreN0I,label=r'$d_R({\cal N}^{(0)},{\cal I})$')
        #plt.semilogx(Frequencies,dFMinAnglestoreN0I,label=r'$d_F({\cal N}^{(0)},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min,AnglestoreN0Ifmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal I})$ from  $d_E( {\cal N}^{(0)},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Icommeasapprxconstsortedmaxdiff_min,AnglestoreN0Icommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal I})$ from  $d_C( {\cal N}^{(0)},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

        # Plot of vairation of |theta| kvec with omgea
        #fig=plt.figure()
        #plt.semilogx(Frequencies,abs(kthetavec[:,0]),label=r'$|\theta k_1|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,1]),label=r'$|\theta k_2|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,2]),label=r'$|\theta k_3|$')

        #plt.xlabel(r'$\omega$ [rad/s]')
        #plt.ylabel(r'$|\theta k_i |$ [rad]')
        #plt.legend()
        #plt.show()

    N0IResults={"Frequencies": Frequencies, "AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max": AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreN0Ifmeasapprxconstsortedmindiff_min": AnglestoreN0Ifmeasapprxconstsortedmindiff_min,"AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min":AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min,\
        "AnglestoreN0Ifmeasapprxconstsortedmindiff_max": AnglestoreN0Ifmeasapprxconstsortedmindiff_max,\
        "AnglestoreN0Icommeasapprxconstsortedmaxdiff_max":AnglestoreN0Icommeasapprxconstsortedmaxdiff_max,"AnglestoreN0Icommeasapprxconstsortedmindiff_min":AnglestoreN0Icommeasapprxconstsortedmindiff_min,\
        "AnglestoreN0Icommeasapprxconstsortedmaxdiff_min":AnglestoreN0Icommeasapprxconstsortedmaxdiff_min,"AnglestoreN0Icommeasapprxconstsortedmindiff_max":AnglestoreN0Icommeasapprxconstsortedmindiff_max,\
        "N0Ifmeasapprx_den_const_max": N0Ifmeasapprx_den_const_max, "N0Ifmeasapprx_den_const_min": N0Ifmeasapprx_den_const_min,\
        "N0Icommeapprx_den_const_max": N0Icommeapprx_den_const_max, "N0Icommeapprx_den_const_min": N0Icommeapprx_den_const_min,\
        "N0Ikthetavec":kthetavec}



########################## N0  and # R

    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    #MinAnglestoreN0R, MaxAnglestoreN0R, dFMinAnglestoreN0R, dFMaxAnglestoreN0R,SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore = MinMaxthetafromQRQI(Frequencies,QN0store,QRstore,UN0store, URstore,MultN0store,MultRstore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultN0store, SortedMultRstore, SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore = SortEigenValues(MultN0store, MultRstore, UN0store, URstore, QN0store, QRstore, Frequencies, sorteigenvalues, N0store, Rstore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreN0Rsortedmaxdiff = AnglesSortedQRQI(SortedQN0store,SortedQRstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Rfmeasfullconstsortedmaxdiff, AnglestoreN0Rfmeasapprxconstsortedmaxdiff_min, AnglestoreN0Rfmeasapprxconstsortedmaxdiff_max, N0Rfmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore, N0store,Rstore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreN0Rcommeasfullconstsortedmaxdiff, AnglestoreN0Rcommeasapprxconstsortedmaxdiff_min, AnglestoreN0Rcommeasapprxconstsortedmaxdiff_max, N0Rcommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore, N0store,Rstore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    SortedMultN0store, SortedMultRstore, SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore = SortEigenValues(MultN0store, MultRstore, UN0store, URstore, QN0store, QRstore, Frequencies, sorteigenvalues, N0store, Rstore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreN0Rsortedmindiff = AnglesSortedQRQI(SortedQN0store,SortedQRstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Rfmeasfullconstsortedmindiff, AnglestoreN0Rfmeasapprxconstsortedmindiff_min,AnglestoreN0Rfmeasapprxconstsortedmindiff_max, N0Rfmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore, N0store,Rstore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    AnglestoreN0Rcommeasfullconstsortedmindiff, AnglestoreN0Rcommeasapprxconstsortedmindiff_min,AnglestoreN0Rcommeasapprxconstsortedmindiff_max, N0Rcommeapprx_den_const_min= Commeasure(sorteigenvalues,SortedUN0store, SortedURstore, SortedQN0store, SortedQRstore, SortedKstore, N0store,Rstore, Frequencies)

    kthetavec=np.zeros((N,3))
    #for n in range(N):
    #    kthetavec[n,0] = -SortedKstore[n,1,2]*MinAnglestoreN0R[n]
    #    kthetavec[n,1] = SortedKstore[n,0,2]*MinAnglestoreN0R[n]
    #    kthetavec[n,2] = -SortedKstore[n,0,1]*MinAnglestoreN0R[n]

    if Figures=="On":
        fig=plt.figure()
        #plt.semilogx(Frequencies,MinAnglestoreN0R,label=r'$d_R({\cal N}^{(0)},{\cal R})$')
        #plt.semilogx(Frequencies,dFMinAnglestoreN0R,label=r'$d_F({\cal N}^{(0)},{\cal R})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Rfmeasapprxconstsortedmaxdiff_min,AnglestoreN0Rfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal R})$ from  $d_E( {\cal N}^{(0)},{\cal R})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Rcommeasapprxconstsortedmaxdiff_min,AnglestoreN0Rcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal R})$ from  $d_C( {\cal N}^{(0)},{\cal R})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

        # Plot of vairation of |theta| kvec with omgea
        #fig=plt.figure()
        #plt.semilogx(Frequencies,abs(kthetavec[:,0]),label=r'$|\theta k_1|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,1]),label=r'$|\theta k_2|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,2]),label=r'$|\theta k_3|$')

        #plt.xlabel(r'$\omega$ [rad/s]')
        #plt.ylabel(r'$|\theta k_i |$ [rad]')
        #plt.legend()
        #plt.show()

    N0RResults={"Frequencies": Frequencies, "AnglestoreN0Rfmeasapprxconstsortedmaxdiff_max": AnglestoreN0Rfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreN0Rfmeasapprxconstsortedmindiff_min": AnglestoreN0Rfmeasapprxconstsortedmindiff_min,"AnglestoreN0Rfmeasapprxconstsortedmaxdiff_min":AnglestoreN0Rfmeasapprxconstsortedmaxdiff_min,\
        "AnglestoreN0Rcommeasapprxconstsortedmaxdiff_max":AnglestoreN0Rcommeasapprxconstsortedmaxdiff_max,
        "AnglestoreN0Rfmeasapprxconstsortedmindiff_max": AnglestoreN0Rfmeasapprxconstsortedmindiff_max,\
        "AnglestoreN0Rcommeasapprxconstsortedmaxdiff_max":AnglestoreN0Rcommeasapprxconstsortedmaxdiff_max,"AnglestoreN0Rcommeasapprxconstsortedmindiff_min":AnglestoreN0Rcommeasapprxconstsortedmindiff_min,\
        "AnglestoreN0Rcommeasapprxconstsortedmaxdiff_min":AnglestoreN0Rcommeasapprxconstsortedmaxdiff_min,"AnglestoreN0Rcommeasapprxconstsortedmindiff_max":AnglestoreN0Rcommeasapprxconstsortedmindiff_max,\
        "N0Rfmeasapprx_den_const_max": N0Rfmeasapprx_den_const_max, "N0Rfmeasapprx_den_const_min": N0Rfmeasapprx_den_const_min,\
        "N0Rcommeapprx_den_const_max": N0Rcommeapprx_den_const_max, "N0Rcommeapprx_den_const_min": N0Rcommeapprx_den_const_min,\
        "N0Rkthetavec":kthetavec}


######################### N0  and # Rtilde

    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    #MinAnglestoreN0Rtilde, MaxAnglestoreN0Rtilde, dFMinAnglestoreN0Rtilde, dFMaxAnglestoreN0Rtilde, SortedUN0store, SortedURtildestore, SortedQN0store, SortedQRtildestore, SortedKstore = MinMaxthetafromQRQI(Frequencies,QN0store,QRtildestore,UN0store, URtildestore,MultN0store,MultRtildestore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    #SortedMultN0store, SortedMultRtildestore, SortedUN0store, SortedURtildestore, SortedQN0store, SortedQRtildestore, SortedKstore = SortEigenValues(MultN0store, MultRtildestore, UN0store, URtildestore, QN0store, QRtildestore, Frequencies, sorteigenvalues, N0store, Rtildestore)
    # Obtain angles for this sorted min-max combination
    AnglestoreN0Rtildesortedmaxdiff = AnglesSortedQRQI(SortedQN0store,SortedQRtildestore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Rtildefmeasfullconstsortedmaxdiff, AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_min,\
    AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_max, N0Rtildefmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedUN0store, SortedURtildestore,
                                                                                                    SortedQN0store, SortedQRtildestore, SortedKstore,
                                                                                                    N0store,Rtildestore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreN0Rtildecommeasfullconstsortedmaxdiff, AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_min, AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_max, N0Rtildecommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedUN0store, SortedURtildestore, SortedQN0store, SortedQRtildestore, SortedKstore, N0store,Rtildestore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    SortedMultN0store, SortedMultRtildestore, SortedUN0store, SortedURtildestore, SortedQN0store, SortedQRtildestore, SortedKstore = SortEigenValues(MultN0store, MultRtildestore, UN0store, URtildestore, QN0store, QRtildestore, Frequencies, sorteigenvalues, N0store, Rtildestore)
    # Obtain angles for this sorted min-max combination
    #AnglestoreN0Rtildesortedmindiff = AnglesSortedQRQI(SortedQN0store,SortedQRtildestore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Rtildefmeasfullconstsortedmindiff, AnglestoreN0Rtildefmeasapprxconstsortedmindiff_min,\
    AnglestoreN0Rtildefmeasapprxconstsortedmindiff_max, N0Rtildefmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedUN0store, SortedURtildestore,
                                                                                                    SortedQN0store, SortedQRtildestore, SortedKstore,
                                                                                                    N0store,Rtildestore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    AnglestoreN0Rtildecommeasfullconstsortedmindiff, AnglestoreN0Rtildecommeasapprxconstsortedmindiff_min,AnglestoreN0Rtildecommeasapprxconstsortedmindiff_max, N0Rtildecommeapprx_den_const_min= Commeasure(sorteigenvalues,SortedUN0store, SortedURtildestore, SortedQN0store, SortedQRtildestore, SortedKstore, N0store,Rtildestore, Frequencies)

    kthetavec=np.zeros((N,3))
    #for n in range(N):
    #    kthetavec[n,0] = -SortedKstore[n,1,2]*MinAnglestoreN0Rtilde[n]
    #    kthetavec[n,1] = SortedKstore[n,0,2]*MinAnglestoreN0Rtilde[n]
    #    kthetavec[n,2] = -SortedKstore[n,0,1]*MinAnglestoreN0Rtilde[n]

    if Figures=="On":
        fig=plt.figure()
        #plt.semilogx(Frequencies,MinAnglestoreN0Rtilde,label=r'$d_R({\cal N}^{(0)},\tilde{\cal R})$')
        #plt.semilogx(Frequencies,dFMinAnglestoreN0Rtilde,label=r'$d_F({\cal N}^{(0)},\tilde{\cal R})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_min,AnglestoreN0Rtildefmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},\tilde{\cal R})$ from  $d_E( {\cal N}^{(0)},\tilde{\cal R})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_min,AnglestoreN0Rtildecommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},\tilde{\cal R})$ from  $d_C( {\cal N}^{(0)},\tilde{\cal R})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

        # Plot of vairation of |theta| kvec with omgea
        #fig=plt.figure()
        #plt.semilogx(Frequencies,abs(kthetavec[:,0]),label=r'$|\theta k_1|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,1]),label=r'$|\theta k_2|$')
        #plt.semilogx(Frequencies,abs(kthetavec[:,2]),label=r'$|\theta k_3|$')

        #plt.xlabel(r'$\omega$ [rad/s]')
        #plt.ylabel(r'$|\theta k_i |$ [rad]')
        #plt.legend()
        #plt.show()

    N0RtildeResults={"Frequencies": Frequencies, "AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_max": AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreN0Rtildefmeasapprxconstsortedmindiff_min": AnglestoreN0Rtildefmeasapprxconstsortedmindiff_min,"AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_min":AnglestoreN0Rtildefmeasapprxconstsortedmaxdiff_min,\
        "AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_max":AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_max,
        "AnglestoreN0Rtildefmeasapprxconstsortedmindiff_max": AnglestoreN0Rtildefmeasapprxconstsortedmindiff_max,\
        "AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_max":AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_max,"AnglestoreN0Rtildecommeasapprxconstsortedmindiff_min":AnglestoreN0Rtildecommeasapprxconstsortedmindiff_min,\
        "AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_min":AnglestoreN0Rtildecommeasapprxconstsortedmaxdiff_min,"AnglestoreN0Rtildecommeasapprxconstsortedmindiff_max":AnglestoreN0Rtildecommeasapprxconstsortedmindiff_max,\
        "N0Rtildefmeasapprx_den_const_max": N0Rtildefmeasapprx_den_const_max, "N0Rtildefmeasapprx_den_const_min": N0Rtildefmeasapprx_den_const_min,\
        "N0Rtildecommeapprx_den_const_max": N0Rtildecommeapprx_den_const_max, "N0Rtildecommeapprx_den_const_min": N0Rtildecommeapprx_den_const_min,\
        "N0Rtildekthetavec":kthetavec}





    return RIResults,RtildeIResults,N0IResults,N0RResults,N0RtildeResults
