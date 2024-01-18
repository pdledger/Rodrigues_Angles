import numpy as np

def EigenValueDecomposition(N0,TensorArray,Frequencies):

    N = len(Frequencies)
    MultRstore=np.zeros(N)
    MultIstore=np.zeros(N)
    MultRtildestore=np.zeros(N)

    URstore=np.zeros((N,3))
    QRstore=np.zeros((N,3,3))
    UIstore=np.zeros((N,3))
    QIstore=np.zeros((N,3,3))
    URtildestore=np.zeros((N,3))
    QRtildestore=np.zeros((N,3,3))

    for n in range(N):
        Mlist = TensorArray[n,:]
        Mten = np.array([[Mlist[0], Mlist[1], Mlist[2]],[Mlist[3], Mlist[4], Mlist[5]],[Mlist[6], Mlist[7], Mlist[8]]])
        Rtilde = np.real(Mten)
        I = np.imag(Mten)
        R = np.real(Mten)-N0

        # Computation of eigenvalues, eigenvector
        uR,VR = np.linalg.eig(R)
        uRtilde,VRtilde = np.linalg.eig(Rtilde)
        uI,VI = np.linalg.eig(I)

    # Possible multiplicities are 1, 2 or 3
    # If they are order in terms of increasing multiplicities we will either have
    #  case 1 : 1 , 1, 1 if they are distinct
    #  case 2 : or 1, 2
    #  case 3 : or 3
        MultR = CheckMult(uR,R)
        MultI = CheckMult(uI,I)
        MultRtilde = CheckMult(uRtilde,Rtilde)

        MultRstore[n] = MultR
        MultIstore[n] = MultI
        MultRtildestore[n] = MultRtilde



        # Store eigenvalues and eigenvectors
        for i in range(3):
            URstore[n,i]=uR[i]
            UIstore[n,i]=uI[i]
            URtildestore[n,i]=uRtilde[i]

            for j in range(3):
                QRstore[n,i,j]=VR[i,j]
                QIstore[n,i,j]=VI[i,j]
                QRtildestore[n,i,j]=VRtilde[i,j]


    return MultRstore, MultIstore, MultRtildestore, URstore, UIstore, URtildestore, QRstore, QIstore, QRtildestore

def CheckMult(u,Tensor):
    # Dynamically adjust tolerance
    Tol=1e-4*np.min(np.abs(u))
    mult=0
    # Determine the multplicity of eigenvalue lambda_i as 3 - rank(R -lambda_i eye(3))
    for i in range(3):
        mult=np.max([mult,3-np.linalg.matrix_rank(Tensor-u[i]*np.eye(3),tol=Tol)])
    return mult
