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
def main(directory,MaxOmega,Figures="On"):
    print("Opening files from this path = ",directory)

    # Read data
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
    FixEvecs="Yes"
    MinAnglestoreRI, MaxAnglestoreRI, dFMinAnglestoreRI, dFMaxAnglestoreRI = MinMaxthetafromQRQI(Frequencies,QRstore,QIstore,URstore, UIstore,MultRstore,MultIstore,FixEvecs)



    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRIsortedmaxdiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRIfmeasfullconstsortedmaxdiff, AnglestoreRIfmeasapprxconstsortedmaxdiff = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)
    AnglestoreRIfmeasfullconstsortedmaxdiff, AnglestoreRIfmeasapprxconstsortedmaxdiff_min,AnglestoreRIfmeasapprxconstsortedmaxdiff_max, RIfmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    #Obtain Com-meauses (approx and exact constant)
    AnglestoreRIcommeasfullconstsortedmaxdiff, AnglestoreRIcommeasapprxconstsortedmaxdiff_min,AnglestoreRIcommeasapprxconstsortedmaxdiff_max, RIcommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRIsortedmindiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRIfmeasfullconstsortedmindiff, AnglestoreRIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)
    AnglestoreRIfmeasfullconstsortedmindiff, AnglestoreRIfmeasapprxconstsortedmindiff_min, AnglestoreRIfmeasapprxconstsortedmindiff_max, RIfmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    #Obtain Com-meauses (approx and exact constant)
    AnglestoreRIcommeasfullconstsortedmindiff, AnglestoreRIcommeasapprxconstsortedmindiff_min,AnglestoreRIcommeasapprxconstsortedmindiff_max, RIcommeapprx_den_const_min = Commeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)


    if Figures=="On":
        fig=plt.figure()
        plt.semilogx(Frequencies,MinAnglestoreRI,label=r'$d_R({\cal R},{\cal I})$')
        plt.semilogx(Frequencies,dFMinAnglestoreRI,label=r'$d_F({\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRIfmeasapprxconstsortedmaxdiff_max,AnglestoreRIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_E({\cal R},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRIcommeasapprxconstsortedmaxdiff_max,AnglestoreRIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_C({\cal R},{\cal I})$ ')

        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

    RIResults={"Frequencies": Frequencies, "MinAnglestoreRI": MinAnglestoreRI, "AnglestoreRIfmeasapprxconstsortedmaxdiff_max": AnglestoreRIfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreRIfmeasapprxconstsortedmindiff_min": AnglestoreRIfmeasapprxconstsortedmindiff_min, "AnglestoreRIfmeasapprxconstsortedmaxdiff_min": AnglestoreRIfmeasapprxconstsortedmaxdiff_min, \
        "AnglestoreRIfmeasapprxconstsortedmindiff_max": AnglestoreRIfmeasapprxconstsortedmindiff_max,
        "dFMinAnglestoreRI":dFMinAnglestoreRI, "dFMaxAnglestoreRI": dFMaxAnglestoreRI,
        "AnglestoreRIcommeasapprxconstsortedmaxdiff_max":AnglestoreRIcommeasapprxconstsortedmaxdiff_max,"AnglestoreRIcommeasapprxconstsortedmindiff_min":AnglestoreRIcommeasapprxconstsortedmindiff_min,\
        "AnglestoreRIcommeasapprxconstsortedmaxdiff_min":AnglestoreRIcommeasapprxconstsortedmaxdiff_min,"AnglestoreRIcommeasapprxconstsortedmindiff_max":AnglestoreRIcommeasapprxconstsortedmindiff_max, \
        "RIfmeasapprx_den_const_max": RIfmeasapprx_den_const_max, "RIfmeasapprx_den_const_min": RIfmeasapprx_den_const_min, \
        "RIcommeapprx_den_const_max":RIcommeapprx_den_const_max, "RIcommeapprx_den_const_min": RIcommeapprx_den_const_min}

    # Plots of minimal and maximal angles and compare with MaxDifference and MinDifference
    if Figures=="On":
        fig=plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(Frequencies,MinAnglestoreRI,label=r'$\theta_{min}$')
        ax1.semilogx(Frequencies,AnglestoreRIsortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta = \theta(Q_R,Q_I)$, $\theta_{max perm} = \theta(Q_R,Q_{I,max perm})$, $\theta_{min perm} = \theta(Q_R,Q_{I,min perm})$ and f measures')
    #ax2 = fig.add_subplot(212)
        ax1.semilogx(Frequencies,MaxAnglestoreRI,label=r'$\theta_{max}$')
        ax1.semilogx(Frequencies,AnglestoreRIsortedmindiff,label=r'$\theta_{min perm}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.legend()
        fig.tight_layout()
        plt.show()

        fig=plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.semilogx(Frequencies,AnglestoreRIsortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')
        ax1.semilogx(Frequencies,np.fmin(AnglestoreRIcommeasapprxconstsortedmaxdiff_max,AnglestoreRIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_C({\cal R},{\cal I})$ ')

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta_{max perm} = \theta(Q_R,Q_{I,max perm})$ and $f$ measures')

        ax2 = fig.add_subplot(212)
        ax2.semilogx(Frequencies,AnglestoreRIsortedmindiff,label=r'$\theta_{min perm}$')
        ax2.semilogx(Frequencies,AnglestoreRIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')
    #ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
        ax2.semilogx(Frequencies,np.fmax(AnglestoreRIcommeasapprxconstsortedmaxdiff_max,AnglestoreRIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_C({\cal R},{\cal I})$ ')

    #ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')
        ax2.set_xlabel(r'$\omega [rad/s]$')
        ax2.set_ylabel(r'$\theta$ [rad]')
        ax2.set_title(r'Comparison of $\theta_{min perm} = \theta(Q_R,Q_{I,min perm})$ and $f$ measures')
        ax1.legend()
        ax2.legend()
        fig.tight_layout()
        plt.show()



    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    MinAnglestoreRtildeI, MaxAnglestoreRtildeI, dFMinAnglestoreRtildeI, dFMaxAnglestoreRtildeI = MinMaxthetafromQRQI(Frequencies,QRtildestore,QIstore,URtildestore, UIstore,MultRtildestore,MultIstore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRtildeIsortedmaxdiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    print('Computing F measure Tilde')
    AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max, RtildeIfmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreRtildeIcommeasfullconstsortedmaxdiff, AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min, AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max, RtildeIcommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRtildeIsortedmindiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max, RtildeIfmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    AnglestoreRtildeIcommeasfullconstsortedmindiff, AnglestoreRtildeIcommeasapprxconstsortedmindiff_min,AnglestoreRtildeIcommeasapprxconstsortedmindiff_max, RtildeIcommeapprx_den_const_min = Commeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)


    if Figures=="On":
        fig=plt.figure()
        plt.semilogx(Frequencies,MinAnglestoreRtildeI,label=r'$d_R(\tilde{\cal R},{\cal I})$')
        plt.semilogx(Frequencies,dFMinAnglestoreRtildeI,label=r'$d_F(\tilde{\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_E(\tilde{\cal R},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max,AnglestoreRtildeIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_C(\tilde{\cal R},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

    RtildeIResults={"Frequencies": Frequencies, "MinAnglestoreRtildeI": MinAnglestoreRtildeI, "AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max": AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreRtildeIfmeasapprxconstsortedmindiff_min": AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,"AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min": AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,
        "AnglestoreRtildeIfmeasapprxconstsortedmindiff_max": AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,
        "AnglestoreRtildeIfmeasfullconstsortedmindiff": AnglestoreRtildeIfmeasfullconstsortedmindiff,\
        "AnglestoreRtildeIfmeasfullconstsortedmaxdiff": AnglestoreRtildeIfmeasfullconstsortedmaxdiff,\
        "dFMinAnglestoreRtildeI":dFMinAnglestoreRtildeI, "dFMaxAnglestoreRtildeI": dFMaxAnglestoreRtildeI,\
        "AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max":AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max,"AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min":AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_min,\
        "AnglestoreRtildeIcommeasapprxconstsortedmindiff_min":AnglestoreRtildeIcommeasapprxconstsortedmindiff_min, "AnglestoreRtildeIcommeasapprxconstsortedmindiff_max":AnglestoreRtildeIcommeasapprxconstsortedmindiff_max, \
        "RtildeIfmeasapprx_den_const_max": RtildeIfmeasapprx_den_const_max, "RtildeIfmeasapprx_den_const_min": RtildeIfmeasapprx_den_const_min,
        "RtildeIcommeapprx_den_const_max":RtildeIcommeapprx_den_const_max, "RtildeIcommeapprx_den_const_min": RtildeIcommeapprx_den_const_min}


    # Plots of minimal and maximal angles and compare with MaxDifference and MinDifference
    if Figures=="On":
        fig=plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(Frequencies,MinAnglestoreRtildeI,label=r'$\theta_{min}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIsortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
    #ax1.set_title(r'Comparison of $\theta_{min} = \min \theta(Q_\tilde{R},Q_I)$  and $\theta_{max perm} = \theta(Q_\tilde{R},Q_{I,max perm})$')
        ax1.semilogx(Frequencies,MaxAnglestoreRtildeI,label=r'$\theta_{max}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIsortedmindiff,label=r'$\theta_{min perm}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax1.set_title(r'Comparison of $\theta = \theta(Q_\tilde{R},Q_I)$, $\theta_{max perm} = \theta(Q_\tilde{R},Q_{I,max perm})$, $\theta_{min perm} = \theta(Q_\tilde{R},Q_{I,min perm})$ and f measures')


    #ax2.set_xlabel(r'$\omega [rad/s]$')
    #ax2.set_ylabel(r'$\theta$ [rad]')
    #ax2.set_title(r'Comparison of $\theta_{max} = \max \theta(Q_\tilde{R},Q_I)$  and $\theta_{min perm} = \theta(Q_\tilde{R},Q_{I,min perm})$')
        ax1.legend()
    #ax2.legend()
        fig.tight_layout()
        plt.show()

        fig=plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.semilogx(Frequencies,AnglestoreRtildeIsortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')
        ax1.semilogx(Frequencies,np.fmin(AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max,AnglestoreRtildeIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_C(\tilde{\cal R},{\cal I})$ ')


        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta_{min} = \theta_{max perm} = \theta(Q_\tilde{R},Q_{I,max perm})$ and $f$ measures')
        ax2 = fig.add_subplot(212)
        ax2.semilogx(Frequencies,AnglestoreRtildeIsortedmindiff,label=r'$\theta_{min perm}$')
        ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')

    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
        ax2.semilogx(Frequencies,np.fmax(AnglestoreRtildeIcommeasapprxconstsortedmaxdiff_max,AnglestoreRtildeIcommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_C(\tilde{\cal R},{\cal I})$ ')

    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax2.set_xlabel(r'$\omega [rad/s]$')
        ax2.set_ylabel(r'$\theta$ [rad]')
        ax2.set_title(r'Comparison of $\theta_{max} = \max \theta(Q_\tilde{R},Q_I)$ and $f$ measures')
        ax1.legend()
        ax2.legend()
        fig.tight_layout()
        plt.show()


    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    MinAnglestoreN0I, MaxAnglestoreN0I, dFMinAnglestoreN0I, dFMaxAnglestoreN0I = MinMaxthetafromQRQI(Frequencies,QN0store,QIstore,UN0store, UIstore,MultN0store,MultIstore,FixEvecs)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultN0store, SortedMultIstore, SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore = SortEigenValues(MultN0store, MultIstore, UN0store, UIstore, QN0store, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreN0Isortedmaxdiff = AnglesSortedQRQI(SortedQN0store,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Ifmeasfullconstsortedmaxdiff, AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min, AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max, N0Ifmeasapprx_den_const_max = Fmeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)

    #Obtain Com meauses (approx and exact constant)
    AnglestoreN0Icommeasfullconstsortedmaxdiff, AnglestoreN0Icommeasapprxconstsortedmaxdiff_min, AnglestoreN0Icommeasapprxconstsortedmaxdiff_max, N0Icommeapprx_den_const_max = Commeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    SortedMultN0store, SortedMultIstore, SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore = SortEigenValues(MultN0store, MultIstore, UN0store, UIstore, QN0store, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreN0Isortedmindiff = AnglesSortedQRQI(SortedQN0store,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    AnglestoreN0Ifmeasfullconstsortedmindiff, AnglestoreN0Ifmeasapprxconstsortedmindiff_min,AnglestoreN0Ifmeasapprxconstsortedmindiff_max, N0Ifmeasapprx_den_const_min = Fmeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)

    #Obtain com meauses (approx and exact constant)
    AnglestoreN0Icommeasfullconstsortedmindiff, AnglestoreN0Icommeasapprxconstsortedmindiff_min,AnglestoreN0Icommeasapprxconstsortedmindiff_max, N0Icommeapprx_den_const_min= Commeasure(sorteigenvalues,SortedUN0store, SortedUIstore, SortedQN0store, SortedQIstore, SortedKstore, N0store,Istore, Frequencies)


    if Figures=="On":
        fig=plt.figure()
        plt.semilogx(Frequencies,MinAnglestoreN0I,label=r'$d_R({\cal N}^{(0)},{\cal I})$')
        plt.semilogx(Frequencies,dFMinAnglestoreN0I,label=r'$d_F({\cal N}^{(0)},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max,AnglestoreN0Ifmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal I})$ from  $d_E( {\cal N}^{(0)},{\cal I})$ ')
        plt.semilogx(Frequencies,np.fmin(AnglestoreN0Icommeasapprxconstsortedmaxdiff_max,AnglestoreN0Icommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R( {\cal N}^{(0)},{\cal I})$ from  $d_C( {\cal N}^{(0)},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

    N0IResults={"Frequencies": Frequencies, "MinAnglestoreN0I": MinAnglestoreN0I, "AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max": AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreN0Ifmeasapprxconstsortedmindiff_min": AnglestoreN0Ifmeasapprxconstsortedmindiff_min,"AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min":AnglestoreN0Ifmeasapprxconstsortedmaxdiff_min,\
        "AnglestoreN0Ifmeasapprxconstsortedmindiff_max": AnglestoreN0Ifmeasapprxconstsortedmindiff_max,\
        "dFMinAnglestoreN0I":dFMinAnglestoreN0I, "dFMaxAnglestoreN0I": dFMaxAnglestoreN0I,\
        "AnglestoreN0Icommeasapprxconstsortedmaxdiff_max":AnglestoreN0Icommeasapprxconstsortedmaxdiff_max,"AnglestoreN0Icommeasapprxconstsortedmindiff_min":AnglestoreN0Icommeasapprxconstsortedmindiff_min,\
        "AnglestoreN0Icommeasapprxconstsortedmaxdiff_min":AnglestoreN0Icommeasapprxconstsortedmaxdiff_min,"AnglestoreN0Icommeasapprxconstsortedmindiff_max":AnglestoreN0Icommeasapprxconstsortedmindiff_max,\
        "N0Ifmeasapprx_den_const_max": N0Ifmeasapprx_den_const_max, "N0Ifmeasapprx_den_const_min": N0Ifmeasapprx_den_const_min,\
        "N0Icommeapprx_den_const_max": N0Icommeapprx_den_const_max, "N0Icommeapprx_den_const_min": N0Icommeapprx_den_const_min}



    # Plots of minimal and maximal angles and compare with MaxDifference and MinDifference
    if Figures=="On":
        fig=plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(Frequencies,MinAnglestoreN0I,label=r'$\theta_{min}$')
        ax1.semilogx(Frequencies,AnglestoreN0Isortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
    #ax1.set_title(r'Comparison of $\theta_{min} = \min \theta(Q_\tilde{R},Q_I)$  and $\theta_{max perm} = \theta(Q_\tilde{R},Q_{I,max perm})$')
        ax1.semilogx(Frequencies,MaxAnglestoreN0I,label=r'$\theta_{max}$')
        ax1.semilogx(Frequencies,AnglestoreN0Isortedmindiff,label=r'$\theta_{min perm}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax1.set_title(r'Comparison of $\theta = \theta(Q_{N},Q_I)$, $\theta_{max perm} = \theta(Q_{N},Q_{I,max perm})$, $\theta_{min perm} = \theta(Q_{N},Q_{I,min perm})$ and f measures')


    #ax2.set_xlabel(r'$\omega [rad/s]$')
    #ax2.set_ylabel(r'$\theta$ [rad]')
    #ax2.set_title(r'Comparison of $\theta_{max} = \max \theta(Q_\tilde{R},Q_I)$  and $\theta_{min perm} = \theta(Q_\tilde{R},Q_{I,min perm})$')
        ax1.legend()
    #ax2.legend()
        fig.tight_layout()
        plt.show()

        fig=plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.semilogx(Frequencies,AnglestoreN0Isortedmaxdiff,label=r'$\theta_{max perm}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasfullconstsortedmaxdiff,label=r'$\theta_{f, max perm, exact C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff,'x',label=r'$\theta_{f, max perm, approx C}$')
    #ax1.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min,'x',label=r'$\theta_{f, max perm, approx C^-}$')
        ax1.semilogx(Frequencies,AnglestoreN0Ifmeasapprxconstsortedmaxdiff_max,'x',label=r'$\theta_{f, max perm, approx C^+}$')
        ax1.semilogx(Frequencies,np.fmin(AnglestoreN0Icommeasapprxconstsortedmaxdiff_max,AnglestoreN0Icommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal N}^{(0)},{\cal I})$ from  $d_C({\cal N}^{(0)},{\cal I})$ ')


        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta_{min} = \theta_{max perm} = \theta(Q_{N},Q_{I,max perm})$ and $f$ measures')
        ax2 = fig.add_subplot(212)
        ax2.semilogx(Frequencies,AnglestoreN0Isortedmindiff,label=r'$\theta_{min perm}$')
        ax2.semilogx(Frequencies,AnglestoreN0Ifmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')

    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax2.semilogx(Frequencies,AnglestoreN0Ifmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
        ax2.semilogx(Frequencies,np.fmax(AnglestoreN0Icommeasapprxconstsortedmaxdiff_max,AnglestoreN0Icommeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal N}^{(0)},{\cal I})$ from  $d_C({\cal N}^{(0)},{\cal I})$ ')

    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax2.set_xlabel(r'$\omega [rad/s]$')
        ax2.set_ylabel(r'$\theta$ [rad]')
        ax2.set_title(r'Comparison of $\theta_{max} = \max \theta(Q_{N},Q_I)$ and $f$ measures')
        ax1.legend()
        ax2.legend()
        fig.tight_layout()
        plt.show()



    return RIResults,RtildeIResults,N0IResults
