import numpy as np
import scipy.linalg
import time

from Rodrigues import *
from StableAngle import *
from CheckOrdering import *

def MinMaxthetafromQRQI(Frequencies,QRstore,QIstore,URstore, UIstore,MultRstore,MultIstore,FixEvecs):

    # For each freuqncy use multiplicties to determine method for
    # obtaining angles
    N = len(Frequencies)
    MinAnglestore = np.zeros(N)
    MaxAnglestore = np.zeros(N)
    dFMinAnglestoreRI = np.zeros(N)
    dFMaxAnglestoreRI = np.zeros(N)


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
        #print(Rmult,Imult)
        if Rmult != Imult :
            print("error different multiplicties for R and I",Rmult,Imult)
            Tol=1e-4*np.min(np.abs(uR))
            #for i in range(3):
            #    print(uR[i],np.linalg.matrix_rank(R-uR[i]*np.eye(3),tol=Tol))
            #Tol=1e-4*np.min(np.abs(uI))
            #for i in range(3):
            #    print(uI[i],np.linalg.matrix_rank(I-uI[i]*np.eye(3),tol=Tol))

            if Rmult==2 or Imult==2:
                Rmult=2
                Imult=2
            else:
                Rmult=1
                Imult=1


        if Rmult==3:
            # This is a case where all 3 eigenvalues are the same. So we know that
            QR=QI
            min_angle, K, tvec = Rodrigues(QR,QI)
            max_angle, K, tvec = Rodrigues(QR,QI)
            dF_min = dFmetric(QR,QI)
            dF_max = dFmetric(QR,QI)

        if Rmult==2:
            # Do a check for the first eigenvector
            angle = np.zeros(6)
            #print(QR,QI,uR,uI)
            #time.sleep(10)
            # Arrange for a mimimal angle
            if FixEvecs=="Yes":
                if np.sign(QR[:,0]@QI[:,0])< 0:
                    QR[:,0] = - QR[:,0]
                QR[:,1:3] = QI[:,1:3]
            else:
                QR, QI, uR, uI = CheckOrdering(QR,QI,uR,uI,0)
            #min_angle=np.arccos(QR[:,0]@QI[:,0])
            #K=np.zeros((3,3))
            #K[1,2]=-1
            #K[2,1]=1
            # Find the combinations of QR and QI that lead to a minimal angle
            #QR, QI, uR = CheckOrdering(QR,QI,uR,0)
            min_angle, K, tvec = Rodrigues(QR,QI)
            #print(min_angle)
            dF_min = dFmetric(QR,QI)
            print(min_angle)
            #time.sleep(10)

            #QR, QI, uR = CheckOrdering(QR,QI,uR,1)
            # Arrange for a maximal angle
            if FixEvecs=="Yes":
                if np.sign(QR[:,0]@QI[:,0])> 0:
                    QR[:,0] = - QR[:,0]
                QR[:,1:3] = QI[:,1:3]
            else:
                QR, QI, uR, uI = CheckOrdering(QR,QI,uR,uI,1)
            max_angle, K, tvec = Rodrigues(QR,QI)
            dF_max = dFmetric(QR,QI)


        if Rmult==1:
            # Find the combinations of QR and QI that lead to a minimal angle
            QR, QI, uR, uI = CheckOrdering(QR,QI,uR,uI,0)
            min_angle, K, tvec = Rodrigues(QR,QI)
            dF_min = dFmetric(QR,QI)
            QR, QI, uR, uI = CheckOrdering(QR,QI,uR,uI,1)
            max_angle, K, tvec = Rodrigues(QR,QI)
            dF_max = dFmetric(QR,QI)

        MinAnglestore[n] = min_angle
        MaxAnglestore[n] = max_angle
        dFMinAnglestoreRI[n] = dF_min
        dFMaxAnglestoreRI[n] = dF_max
        if n > 5 and n < 8:
            print(n,Rmult)
            print(uR,uI)
            print(QR,QI)
            print(min_angle,max_angle)
            print(dF_min,dF_max)





    return MinAnglestore, MaxAnglestore, dFMinAnglestoreRI, dFMaxAnglestoreRI

def dFmetric(QR,QI):
    return np.sqrt(np.trace((QR-QI)@np.transpose(QR-QI)))
