import numpy as np
import time
def Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies):
    N=len(Frequencies)
    Fexactconst = np.zeros(N,dtype=np.longdouble)
    Fexactconst3 = np.zeros(N,dtype=np.longdouble)
    Fexactconst4 = np.zeros(N,dtype=np.longdouble)


    #Fapproxconst = np.zeros(N)
    Fapproxconst_min = np.zeros(N,dtype=np.longdouble)
    Fapproxconst_max = np.zeros(N,dtype=np.longdouble)
    den_const=np.zeros(N)
    # First compute F using the exact constant
    for n in range(N):
        QR = np.zeros((3,3),dtype=np.longdouble )
        QI = np.zeros((3,3),dtype=np.longdouble)
        R=np.zeros((3,3),dtype=np.longdouble)
        I=np.zeros((3,3),dtype=np.longdouble)
        uR = np.zeros(3,dtype=np.longdouble)
        uI = np.zeros(3,dtype=np.longdouble)
        K = np.zeros((3,3),dtype=np.longdouble)
        for i in range(3):
            uR[i] = SortedURstore[n,i]
            uI[i] = SortedUIstore[n,i]
            for j in range(3):
                QR[i,j] = SortedQRstore[n,i,j]
                QI[i,j] = SortedQIstore[n,i,j]
                K[i,j] = SortedKstore[n,i,j]
                R[i,j] = Rstore[n,i,j]
                I[i,j] = Istore[n,i,j]
        diffeig=0.
        for i in range(3):
            diffeig+=(uR[i]-uI[i])**2
        normalisation =np.abs( np.trace(K@K@np.diag((uR))@np.diag((uI))))- np.abs(np.trace(K@np.diag((uR))@K@np.diag((uI))))
        Fexactconst[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation)
        #normalisation =np.abs( np.trace(K@K@np.diag((uR))@K@np.diag((uI))))- np.abs(np.trace(K@K@np.diag((uI))@K@np.diag((uR))))
        #Fexactconst3[n] = 2*np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation)
        #normalisation =np.abs( np.trace(K@K@np.diag((uR))@K@K@np.diag((uI))))
        #Fexactconst4[n] = 4*np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation)

    Fexactconst= np.sqrt(Fexactconst)
    #Fexactconst3= (Fexactconst3)**(1/3)
    #Fexactconst4= (Fexactconst4)**(1/4)




    # Next compute F using the apprx constant without K, QR or QI
    for n in range(N):
        R=np.zeros((3,3),dtype=np.longdouble)
        I=np.zeros((3,3),dtype=np.longdouble)
        uR = np.zeros(3,dtype=np.longdouble)
        uI = np.zeros(3,dtype=np.longdouble)
        for i in range(3):
            uR[i] = SortedURstore[n,i]
            uI[i] = SortedUIstore[n,i]
            for j in range(3):
                R[i,j] = Rstore[n,i,j]
                I[i,j] = Istore[n,i,j]
        diffeig=0.
        #print(Frequencies[n],uR,uI)
        for i in range(3):
            diffeig=diffeig+(uR[i]-uI[i])**2

        #normalisationapprox=0
        #for j in range(3):
        #    normalisationapprox-=2*(uI[j]*uR[j])/np.sqrt(3)
        #normalisationapprox+=((uI[1]*uR[0])+(uI[2]*uR[0])+(uI[0]*uR[1])+(uI[2]*uR[1])+(uI[0]*uR[2])+(uI[1]*uR[2]))/np.sqrt(3)
        evlist = np.zeros(3,dtype=np.longdouble)
        #evlist[0]= -uI[1]*uR[1]-uI[2]*uR[2]+uI[2]*uR[1]+uI[1]*uR[2]
        #evlist[1]= -uI[0]*uR[0]-uI[2]*uR[2]+uI[2]*uR[0]+uI[0]*uR[2]
        #evlist[2]= -uI[1]*uR[1]-uI[0]*uR[0]+uI[1]*uR[0]+uI[0]*uR[1]
        evlist[0]= - (uI[1]-uI[2])*(uR[1]-uR[2])
        evlist[1]= - (uI[0]-uI[2])*(uR[0]-uR[2])
        evlist[2]= - (uI[0]-uI[1])*(uR[0]-uR[1])


        normalisation_min = np.min(np.abs(evlist))
        normalisation_max = np.max(np.abs(evlist))
        Tol=1e-6
        #Fapproxconst_min[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_min)
        #Fapproxconst_max[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_max)
        Fapproxconst_min[n] = np.min([(np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig)) / (normalisation_min),
        (np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig)) / (normalisation_max)])
        Fapproxconst_max[n] = np.max([(np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig)) / (normalisation_min),
        (np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig)) / (normalisation_max)])

        #Fapproxconst_min[n] = np.min([np.abs(np.linalg.norm(R/ np.sqrt(normalisation_min) - I /np.sqrt(normalisation_min) ,ord='fro')**2 - diffeig / normalisation_min),
        #np.abs(np.linalg.norm(R/ np.sqrt(normalisation_max)- I/ np.sqrt(normalisation_max) ,ord='fro')**2 - diffeig / normalisation_max)])
        #Fapproxconst_max[n] = np.max([np.abs(np.linalg.norm(R/ np.sqrt(normalisation_min)-I/ np.sqrt(normalisation_min),ord='fro')**2 - diffeig / normalisation_min),
        #np.abs(np.linalg.norm(R/ np.sqrt(normalisation_max) - I/ np.sqrt(normalisation_max),ord='fro')**2 - diffeig / normalisation_max)])


        #print(Frequencies[n],(np.linalg.norm(R-I,ord='fro')**2 - diffeig),normalisation_min,normalisation_max)
        #if Fapproxconst_min[n] <0:
        #    Fapproxconst_min[n] = - Fapproxconst_min[n]
        #if Fapproxconst_max[n] <0:
        #    Fapproxconst_max[n] = - Fapproxconst_max[n]
        #print("fmes",Frequencies[n],(np.linalg.norm(R-I,ord='fro')**2 - diffeig),normalisation_min,normalisation_max)
        #print(R-np.transpose(R),I-np.transpose(I))
        #print(R,I)
        #print(uR,uI)
        #uR,VR = np.linalg.eig(R.astype(dtype=float))
        #uI,VI = np.linalg.eig(I.astype(dtype=float))
        #print(uR,uI)
        #time.sleep(10)

        # if abs(normalisation_min/(np.linalg.norm(uI)*np.linalg.norm(uR))) > Tol**2:
        #     Fapproxconst_min[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_min)
        # else:
        #     print(normalisation_min/(np.linalg.norm(uI)*np.linalg.norm(uR)))
        #     Fapproxconst_min[n] = 0
        # Tol=1e-6
        # if abs(normalisation_max/(np.linalg.norm(uI)*np.linalg.norm(uR))) > Tol**2:
        #     Fapproxconst_max[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_max)
        # else:
        #     print(normalisation_max/(np.linalg.norm(uI)*np.linalg.norm(uR)))
        #     Fapproxconst_max[n] = 0
        den_const[n]=np.min([np.sqrt(np.abs(normalisation_min)),np.sqrt(np.abs(normalisation_max))])
        #print(Frequencies[n],Fapproxconst_min[n],np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig),np.abs(normalisation_min))
        #print(Frequencies[n],Fapproxconst_max[n],np.abs(normalisation_max))
#
#        Fapproxconst[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisationapprox)
#    Fapproxconst= np.sqrt(Fapproxconst)
    Fapproxconst_min=np.sqrt(Fapproxconst_min)
    Fapproxconst_max=np.sqrt(Fapproxconst_max)

    print(Fapproxconst_min)

    return Fexactconst,Fapproxconst_min,Fapproxconst_max,den_const
    #Fapproxconst
