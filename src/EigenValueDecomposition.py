import numpy as np

def EigenValueDecomposition(N0,TensorArray,Frequencies):

    N = len(Frequencies)
    MultRstore=np.zeros(N, dtype=np.longdouble)
    MultIstore=np.zeros(N, dtype=np.longdouble)
    MultRtildestore=np.zeros(N, dtype=np.longdouble)
    MultN0store=np.zeros(N, dtype=np.longdouble) # repeated storage of N0 multplicities, eigenvectors, eigenvalues etc

    URstore=np.zeros((N,3), dtype=np.longdouble)
    QRstore=np.zeros((N,3,3), dtype=np.longdouble)
    UIstore=np.zeros((N,3), dtype=np.longdouble)
    QIstore=np.zeros((N,3,3), dtype=np.longdouble)
    UN0store=np.zeros((N,3), dtype=np.longdouble)
    QN0store=np.zeros((N,3,3), dtype=np.longdouble)
    URtildestore=np.zeros((N,3), dtype=np.longdouble)
    QRtildestore=np.zeros((N,3,3), dtype=np.longdouble)

    for n in range(N):
        Mlist = TensorArray[n,:]
        Mten = np.array([[Mlist[0], Mlist[1], Mlist[2]],[Mlist[3], Mlist[4], Mlist[5]],[Mlist[6], Mlist[7], Mlist[8]]], dtype=np.longcomplex)
        Rtilde = np.real(Mten)
        I = np.imag(Mten)
        R = (np.real(Mten)-N0)

        # Computation of eigenvalues, eigenvector
        uR,VR = np.linalg.eig(R.astype(dtype=float))
        uRtilde,VRtilde = np.linalg.eig(Rtilde.astype(dtype=float))
        uI,VI = np.linalg.eig(I.astype(dtype=float))
        uN0,VN0 = np.linalg.eig(N0.astype(dtype=float))


    # Possible multiplicities are 1, 2 or 3
    # If they are order in terms of increasing multiplicities we will either have
    #  case 1 : 1 , 1, 1 if they are distinct
    #  case 2 : or 1, 2
    #  case 3 : or 3
        MultR = CheckMult(uR,R)
        MultI = CheckMult(uI,I)
        MultRtilde = CheckMult(uRtilde,Rtilde)
        MultN0 = CheckMult(uN0,N0)


        if n > 5 and n < 8:
            print(n)
            print("R",uR,VR,MultR)
            print("Rtilde",uRtilde,VRtilde,MultRtilde)
            print("I",uI,VI,MultI)


        MultRstore[n] = MultR
        MultIstore[n] = MultI
        MultRtildestore[n] = MultRtilde
        MultN0store[n] = MultN0



        # Store eigenvalues and eigenvectors
        for i in range(3):
            URstore[n,i]=uR[i]
            UIstore[n,i]=uI[i]
            URtildestore[n,i]=uRtilde[i]
            UN0store[n,i]=uN0[i]

            for j in range(3):
                QRstore[n,i,j]=VR[i,j]
                QIstore[n,i,j]=VI[i,j]
                QRtildestore[n,i,j]=VRtilde[i,j]
                QN0store[n,i,j]=VN0[i,j]

    return MultRstore, MultIstore, MultRtildestore, MultN0store, URstore, UIstore, URtildestore, UN0store, QRstore, QIstore, QRtildestore, QN0store

def CheckMult(u,Tensor):
    # Dynamically adjust tolerance
#    Tol=1e-4*np.min(np.abs(u))
#    mult=0
#    # Determine the multplicity of eigenvalue lambda_i as 3 - rank(R -lambda_i eye(3))
#    for i in range(3):
#        mult=np.max([mult,3-np.linalg.matrix_rank(Tensor-u[i]*np.eye(3),tol=Tol)])
#
    Tol=5e-4
    mult=1
    for i in range(3):
        for j in range(i+1,3):
            if i != j:
                if np.abs(u[i]-u[j])/abs(u[i]) < Tol:
                    mult=mult+1
                    if mult > 3:
                        mult =3

    return mult
