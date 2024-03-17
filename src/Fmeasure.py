import numpy as np
def Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies):
    N=len(Frequencies)
    Fexactconst = np.zeros(N)
    #Fapproxconst = np.zeros(N)
    Fapproxconst_min = np.zeros(N)
    Fapproxconst_max = np.zeros(N)
    # First compute F using the exact constant
    for n in range(N):
        QR = np.zeros((3,3))
        QI = np.zeros((3,3))
        R=np.zeros((3,3))
        I=np.zeros((3,3))
        uR = np.zeros(3)
        uI = np.zeros(3)
        K = np.zeros((3,3))
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
            diffeig+=np.abs(uR[i]-uI[i])**2
        normalisation =np.abs( np.trace(K@K@np.diag((uR))@np.diag((uI))))- np.abs(np.trace(K@np.diag((uR))@K@np.diag((uI))))
        Fexactconst[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation)
    Fexactconst= np.sqrt(Fexactconst)


    # Next compute F using the apprx constant without K, QR or QI
    for n in range(N):
        R=np.zeros((3,3))
        I=np.zeros((3,3))
        uR = np.zeros(3)
        uI = np.zeros(3)
        for i in range(3):
            uR[i] = SortedURstore[n,i]
            uI[i] = SortedUIstore[n,i]
            for j in range(3):
                R[i,j] = Rstore[n,i,j]
                I[i,j] = Istore[n,i,j]
        diffeig=0.
        for i in range(3):
            diffeig+=np.abs(uR[i]-uI[i])**2

        #normalisationapprox=0
        #for j in range(3):
        #    normalisationapprox-=2*(uI[j]*uR[j])/np.sqrt(3)
        #normalisationapprox+=((uI[1]*uR[0])+(uI[2]*uR[0])+(uI[0]*uR[1])+(uI[2]*uR[1])+(uI[0]*uR[2])+(uI[1]*uR[2]))/np.sqrt(3)
        evlist = np.zeros(3)
        #evlist[0]= -uI[1]*uR[1]-uI[2]*uR[2]+uI[2]*uR[1]+uI[1]*uR[2]
        #evlist[1]= -uI[0]*uR[0]-uI[2]*uR[2]+uI[2]*uR[0]+uI[0]*uR[2]
        #evlist[2]= -uI[1]*uR[1]-uI[0]*uR[0]+uI[1]*uR[0]+uI[0]*uR[1]
        evlist[0]= - (uI[1]-uI[2])*(uR[1]-uR[2])
        evlist[1]= - (uI[0]-uI[2])*(uR[0]-uR[2])
        evlist[2]= - (uI[0]-uI[1])*(uR[0]-uR[1])


        normalisation_min = np.min(evlist)
        normalisation_max = np.max(evlist)
        Tol=1e-6
        if abs(normalisation_min/(np.linalg.norm(uI)*np.linalg.norm(uR))) > Tol**2:
            Fapproxconst_min[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_min)
        else:
            print(normalisation_min/(np.linalg.norm(uI)*np.linalg.norm(uR)))
            Fapproxconst_min[n] = 0
        Tol=1e-6
        if abs(normalisation_max/(np.linalg.norm(uI)*np.linalg.norm(uR))) > Tol**2:
            Fapproxconst_max[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisation_max)
        else:
            print(normalisation_max/(np.linalg.norm(uI)*np.linalg.norm(uR)))
            Fapproxconst_max[n] = 0

        #print(Frequencies[n],Fapproxconst_min[n],np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig),np.abs(normalisation_min))
        #print(Frequencies[n],Fapproxconst_max[n],np.abs(normalisation_max))
#
#        Fapproxconst[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisationapprox)
#    Fapproxconst= np.sqrt(Fapproxconst)
    Fapproxconst_min=np.sqrt(Fapproxconst_min)
    Fapproxconst_max=np.sqrt(Fapproxconst_max)

    return Fexactconst,Fapproxconst_min,Fapproxconst_max
    #Fapproxconst
