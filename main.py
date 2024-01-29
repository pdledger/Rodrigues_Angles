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
    for n in range(N):
        if Frequencies[n] > MaxOmega:
            Ntarget=n
            break

    # allow for frequency array less than max omega.
    if MaxOmega > np.max(Frequencies):
        Ntarget = -1

    TensorArray = TensorArray[:Ntarget,:]
    Frequencies = Frequencies[:Ntarget]
    N=Ntarget

    Rstore,Istore,Rtildestore = SplitTensor(TensorArray,Frequencies,N0)

    # Determine eigenvalue decompositions of N0, R, I, Rtilde (no sorting applied), and their multiplicities
    MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore = EigenValueDecomposition(N0,TensorArray,Frequencies)

    # Determine the maximal and minimal angles from QR and QI also output d_F metric for these orderings
    MinAnglestoreRI, MaxAnglestoreRI, dFMinAnglestoreRI, dFMaxAnglestoreRI = MinMaxthetafromQRQI(Frequencies,QRstore,QIstore,URstore, UIstore,MultRstore,MultIstore)



    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRIsortedmaxdiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRIfmeasfullconstsortedmaxdiff, AnglestoreRIfmeasapprxconstsortedmaxdiff = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)
    AnglestoreRIfmeasfullconstsortedmaxdiff, AnglestoreRIfmeasapprxconstsortedmaxdiff_min,AnglestoreRIfmeasapprxconstsortedmaxdiff_max = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore = SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRIsortedmindiff = AnglesSortedQRQI(SortedQRstore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRIfmeasfullconstsortedmindiff, AnglestoreRIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)
    AnglestoreRIfmeasfullconstsortedmindiff, AnglestoreRIfmeasapprxconstsortedmindiff_min, AnglestoreRIfmeasapprxconstsortedmindiff_max = Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies)

    if Figures=="On":
        fig=plt.figure()
        plt.semilogx(Frequencies,MinAnglestoreRI,label=r'$d_R({\cal R},{\cal I})$')
        plt.semilogx(Frequencies,dFMinAnglestoreRI,label=r'$d_F({\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRIfmeasapprxconstsortedmaxdiff_max,AnglestoreRIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R({\cal R},{\cal I})$ from  $d_E({\cal R},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

    RIResults={"Frequencies": Frequencies, "MinAnglestoreRI": MinAnglestoreRI, "AnglestoreRIfmeasapprxconstsortedmaxdiff_max": AnglestoreRIfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreRIfmeasapprxconstsortedmindiff_min": AnglestoreRIfmeasapprxconstsortedmindiff_min,"dFMinAnglestoreRI":dFMinAnglestoreRI, "dFMaxAnglestoreRI": dFMaxAnglestoreRI}

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
        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta_{max perm} = \theta(Q_R,Q_{I,max perm})$ and $f$ measures')

        ax2 = fig.add_subplot(212)
        ax2.semilogx(Frequencies,AnglestoreRIsortedmindiff,label=r'$\theta_{min perm}$')
        ax2.semilogx(Frequencies,AnglestoreRIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')
    #ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
    #ax2.semilogx(Frequencies,AnglestoreRIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')
        ax2.set_xlabel(r'$\omega [rad/s]$')
        ax2.set_ylabel(r'$\theta$ [rad]')
        ax2.set_title(r'Comparison of $\theta_{min perm} = \theta(Q_R,Q_{I,min perm})$ and $f$ measures')
        ax1.legend()
        ax2.legend()
        fig.tight_layout()
        plt.show()



    # Determine the maximal and minimal angles from QRtilde and QI also output d_F metric for these orderings
    MinAnglestoreRtildeI, MaxAnglestoreRtildeI, dFMinAnglestoreRtildeI, dFMaxAnglestoreRtildeI = MinMaxthetafromQRQI(Frequencies,QRtildestore,QIstore,URtildestore, UIstore,MultRtildestore,MultIstore)

    # Sort eigenvalues (and eigenvectors) so that || Lambda_Rtilde - Lambda_I || is maximal
    sorteigenvalues="MaxDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRtildeIsortedmaxdiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    AnglestoreRtildeIfmeasfullconstsortedmaxdiff, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_min, AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)


    # Sort eigenvalues (and eigenvectors) so that || Lambda_R - Lambda_I || is minimal
    sorteigenvalues="MinDifference"
    #SortedMultRstore, SortedMultIstore, SortedMultRtildestore, SortedURstore, SortedUIstore, SortedURtildestore, SortedQRstore, SortedQIstore, SortedQRtildestore = SortEigenValues(MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore, Frequencies)
    SortedMultRtildestore, SortedMultIstore, SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore = SortEigenValues(MultRtildestore, MultIstore, URtildestore, UIstore, QRtildestore, QIstore, Frequencies, sorteigenvalues)
    # Obtain angles for this sorted min-max combination
    AnglestoreRtildeIsortedmindiff = AnglesSortedQRQI(SortedQRtildestore,SortedQIstore,Frequencies)

    #Obtain f meauses (approx and exact constant)
    #AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)
    AnglestoreRtildeIfmeasfullconstsortedmindiff, AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max = Fmeasure(sorteigenvalues,SortedURtildestore, SortedUIstore, SortedQRtildestore, SortedQIstore, SortedKstore, Rtildestore,Istore, Frequencies)

    if Figures=="On":
        fig=plt.figure()
        plt.semilogx(Frequencies,MinAnglestoreRtildeI,label=r'$d_R(\tilde{\cal R},{\cal I})$')
        plt.semilogx(Frequencies,dFMinAnglestoreRtildeI,label=r'$d_F(\tilde{\cal R},{\cal I})$')
        plt.semilogx(Frequencies,np.fmin(AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min),'x',label=r'Approx $d_R(\tilde{\cal R},{\cal I})$ from  $d_E(\tilde{\cal R},{\cal I})$ ')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\theta$ [rad]')
        plt.legend()
        plt.show()

    RtildeIResults={"Frequencies": Frequencies, "MinAnglestoreRtildeI": MinAnglestoreRtildeI, "AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max": AnglestoreRtildeIfmeasapprxconstsortedmaxdiff_max, \
        "AnglestoreRtildeIfmeasapprxconstsortedmindiff_min": AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,"dFMinAnglestoreRtildeI":dFMinAnglestoreRtildeI, "dFMaxAnglestoreRtildeI": dFMaxAnglestoreRtildeI}


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

        ax1.set_xlabel(r'$\omega [rad/s]$')
        ax1.set_ylabel(r'$\theta$ [rad]')
        ax1.set_title(r'Comparison of $\theta_{min} = \theta_{max perm} = \theta(Q_\tilde{R},Q_{I,max perm})$ and $f$ measures')
        ax2 = fig.add_subplot(212)
        ax2.semilogx(Frequencies,AnglestoreRtildeIsortedmindiff,label=r'$\theta_{min perm}$')
        ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasfullconstsortedmindiff,label=r'$\theta_{f, min perm, exact C}$')

    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff,'x',label=r'$\theta_{f, min perm, approx C}$')
        ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_min,'x',label=r'$\theta_{f, min perm, approx C^-}$')
    #ax2.semilogx(Frequencies,AnglestoreRtildeIfmeasapprxconstsortedmindiff_max,'x',label=r'$\theta_{f, min perm, approx C^+}$')

        ax2.set_xlabel(r'$\omega [rad/s]$')
        ax2.set_ylabel(r'$\theta$ [rad]')
        ax2.set_title(r'Comparison of $\theta_{max} = \max \theta(Q_\tilde{R},Q_I)$ and $f$ measures')
        ax1.legend()
        ax2.legend()
        fig.tight_layout()
        plt.show()



    return RIResults,RtildeIResults
