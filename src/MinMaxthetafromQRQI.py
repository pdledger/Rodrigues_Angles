import numpy as np
import scipy.linalg

from Rodrigues import *
from StableAngle import *
from CheckOrdering import *

def MinMaxthetafromQRQI(Frequencies,QRstore,QIstore,URstore, UIstore,MultRstore,MultIstore):

    # For each freuqncy use multiplicties to determine method for
    # obtaining angles
    N = len(Frequencies)
    MinAnglestore = np.zeros(N)
    MaxAnglestore = np.zeros(N)


    for n in range(N):
        QR = np.zeros((3,3))
        QI = np.zeros((3,3))
        uR =np.zeros(3)
        for i in range(3):
            uR[i]=URstore[n,i]
            for j in range(3):
                QR[i,j] = QRstore[n,i,j]
                QI[i,j] = QIstore[n,i,j]
        Rmult = MultRstore[n]
        Imult = MultIstore[n]

        if Rmult != Imult :
            print("error different multiplicties for R and I",Rmult,Imult)
            Tol=1e-4*np.min(np.abs(uR))
            for i in range(3):
                print(uR[i],np.linalg.matrix_rank(R-uR[i]*np.eye(3),tol=Tol))
            Tol=1e-4*np.min(np.abs(uI))
            for i in range(3):
                print(uI[i],np.linalg.matrix_rank(I-uI[i]*np.eye(3),tol=Tol))

            if Rmult==2 or Imult==2:
                Rmult=2
                Imult=2
            else:
                Rmult=1
                Imult=1


        if Rmult==3:
            # This is a case where all 3 eigenvalues are the same. So we know that
            QR=QI

        if Rmult==2:
            # Do a check for the first eigenvector
            angle = np.zeros(6)
            # Arrange for a mimimal angle
            if np.sign(QR[:,0]@QI[:,0])< 0:
                QR[:,0] = - QR[:,0]
            QR[:,1:3] = QI[:,1:3]

            min_angle, K, tvec = Rodrigues(QR,QI)

            # Arrange for a maximal angle
            if np.sign(QR[:,0]@QI[:,0])> 0:
                QR[:,0] = - QR[:,0]
            QR[:,1:3] = QI[:,1:3]

            max_angle, K, tvec = Rodrigues(QR,QI)


        if Rmult==1:
            # Find the combinations of QR and QI that lead to a minimal angle
            QR, QI, uR = CheckOrdering(QR,QI,uR,0)
            min_angle, K, tvec = Rodrigues(QR,QI)
            QR, QI, uR = CheckOrdering(QR,QI,uR,1)
            max_angle, K, tvec = Rodrigues(QR,QI)

        MinAnglestore[n] = min_angle
        MaxAnglestore[n] = max_angle




    return MinAnglestore, MaxAnglestore
