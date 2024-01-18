import numpy as np
def Fmeasure(sorteigenvalues,SortedURstore, SortedUIstore, SortedQRstore, SortedQIstore, SortedKstore, Rstore,Istore, Frequencies):
    N=len(Frequencies)
    Fexactconst = np.zeros(N)
    Fapproxconst = np.zeros(N)

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

        normalisationapprox=0
        for j in range(3):
            normalisationapprox-=2*(uI[j]*uR[j])
        normalisationapprox+=(uI[1]*uR[0])+(uI[2]*uR[0])+(uI[0]*uR[1])+(uI[2]*uR[1])+(uI[0]*uR[2])+(uI[1]*uR[2])

        Fapproxconst[n] = np.abs(np.linalg.norm(R-I,ord='fro')**2 - diffeig) / np.abs(normalisationapprox)
    Fapproxconst= np.sqrt(Fapproxconst)
    return Fexactconst,Fapproxconst
