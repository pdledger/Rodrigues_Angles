import numpy as np

from Rodrigues import *
def SortEigenValues(MultRstore, MultIstore, URstore, UIstore, QRstore, QIstore, Frequencies, sorteigenvalues):
    N = len(Frequencies)
    # Prepare sorted values by making a copy (note multiplicties don't change)
    SortedMultRstore=np.copy(MultRstore)
    SortedMultIstore=np.copy(MultIstore)
    SortedURstore =np.copy(URstore)
    SortedUIstore =np.copy(UIstore)
    SortedQRstore = np.copy(QRstore)
    SortedQIstore = np.copy(QIstore)
    SortedKstore = np.zeros((N,3,3))

    Perm = np.array([[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]])
    sign=np.array([[1,1,1], \
                        [-1,1,1], \
                        [-1,-1,1], \
                        [-1,-1,-1], \
                        [-1,1,-1], \
                        [1,-1,-1], \
                        [1,-1,1], \
                        [1,1,-1]])


    for n in range(N):
        QR = np.zeros((3,3))
        QI = np.zeros((3,3))
        uR =np.zeros(3)
        uI = np.zeros(3)
        for i in range(3):
            uR[i]=URstore[n,i]
            uI[i]=UIstore[n,i]

            for j in range(3):
                QR[i,j] = QRstore[n,i,j]
                QI[i,j] = QIstore[n,i,j]
        Rmult = MultRstore[n]
        Imult = MultIstore[n]

        if sorteigenvalues=="MinDifference":
            # Find min combination
            diffeig=1e10
        elif sorteigenvalues=="MaxDifference":
            # Find max combination
            diffeig=0.
        for m in range(6):
            sum=0.
            ind=Perm[m,:]
            for i in range(3):
                sum = sum+ abs(uR[i]-uI[ind[i]-1])**2
            check = False
            if sorteigenvalues=="MinDifference" and sum < diffeig:
                check = True
            elif sorteigenvalues=="MaxDifference" and sum > diffeig:
                check = True
            if check==True:
                diffeig = sum
                puI=np.zeros(3)
                #S=np.zeros((3,3))
                for i in range(3):
                    puI[i]=uI[ind[i]-1]

                thetaopt=1e10
                Kopt=np.zeros((3,3))
                for k in range(8):
                    QIordsign=np.zeros((3,3))
                    for j in range(3):
                        QIordsign[:,j] = sign[k,j]*QI[:,ind[j]-1]
                    if np.linalg.det(np.transpose(QR)@QIordsign)> 0:
                        # Only do for valid rotation matrices with det(R) + ve (=1)
                        theta, K, Tvec= Rodrigues(QR, QIordsign)
                        if theta < thetaopt:
                            thetaopt=theta
                            Kopt=np.copy(K)
                            uRopt=np.copy(uR)
                            uIopt =np.copy(puI)
                            QRopt = np.copy(QR)
                            QIopt = np.copy(QIordsign)
        # Store the optimal combination
        for i in range(3):
            SortedURstore[n,i]=uRopt[i]
            SortedUIstore[n,i]=uIopt[i]

            for j in range(3):
                SortedQRstore[n,i,j] = QRopt[i,j]
                SortedQIstore[n,i,j] = QIopt[i,j]
                SortedKstore[n,i,j] = Kopt[i,j]


    return SortedMultRstore, SortedMultIstore, SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore
