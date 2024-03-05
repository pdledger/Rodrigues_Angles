import numpy as np
import time

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
        # nb output from np.linalg.eig is not necessarily in order


    # Possible multiplicities are 1, 2 or 3
    # If they are order in terms of increasing multiplicities we will either have
    #  case 1 : 1 , 1, 1 if they are distinct
    #  case 2 : or 1, 2
    #  case 3 : or 3
        MultR = CheckMult(uR,R)
        MultI = CheckMult(uI,I)
        MultRtilde = CheckMult(uRtilde,Rtilde)

        # Reorder eigenvalues/eigenvectors depending on multiplicity
        uR,VR,MultR=Reorder(uR,VR,MultR)
        uI,VI,MultI=Reorder(uI,VI,MultI)
        uRtilde,VRtilde,MultRtilde=Reorder(uRtilde,VRtilde,MultRtilde)


        MultRstore[n] = np.max(MultR)
        MultIstore[n] = np.max(MultI)
        MultRtildestore[n] = np.max(MultRtilde)

        # Check for reliablity of computed Eigenvalues
        VR=checkreliable(uR,"R",VR)
        VI=checkreliable(uI,"I",VI)
        VRtilde=checkreliable(uRtilde,"Rtilde",VRtilde)



        #print("Eigenvalues,Eignvectors")
        #if n<4:
        #    print(uR,uRtilde,uI)
        #    print(VR)
        #    print(VRtilde)
        #    print(VI)
        #    print(np.max(MultR),np.max(MultRtilde), np.max(MultI))
        #    time.sleep(1)


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
    #print(u)
    mult=np.zeros(3)
    for i in range(3):
        #print(3-np.linalg.matrix_rank(Tensor-u[i]*np.eye(3),tol=Tol))
        mult[i]=3-np.linalg.matrix_rank(Tensor-u[i]*np.eye(3),tol=Tol)
    return mult

def Reorder(u,V,Mult):
    # Choose to re-order eigenvalues/eigenvectors so that they are ordered in increasing multplicity.
    # Possibilities 1)
    # All distinct and have multplicity 1
    # Two eigenvalues have multplicity 2 (swap maybe required.)
    # All eigenvalues the same and have multplicity 3
    Maxmult=np.max(Mult)
    if Maxmult==2:
        # Set the eigenvalue/eigenvector-pair to have multiplicity 1
        newcol0=np.argmin(Mult)
        if newcol0 !=0:
            oldev=np.copy(u[0])
            u[0]=np.copy(u[newcol0])
            u[newcol0]=oldev
            oldevec=np.copy(V[0,:])
            V[0,:]=np.copy(V[newcol0,:])
            V[newcol0,:]=oldevec
            oldmult=np.copy(Mult[0])
            Mult[0]=np.copy(Mult[newcol0])
            Mult[newcol0] =oldmult
            #print(u)
            #print(V)
            #print(Mult)
    return u,V,Mult

def checkreliable(u,Tensor,V):
    Perm = np.array([[0,1],[0,2],[1,2]])
    Tol=1e-2
    for n in range(3):
        p=Perm[n,:]
        diff=np.abs(u[p[0]]-u[p[1]])/np.sqrt(u[p[0]]**2+u[p[1]]**2)
        if diff < Tol:
            print("Eigenvalues closely spaced for",Tensor,diff,Tol)
            print("Permutation",n)
            V=np.eye(3)
        #print(diff,Tensor)
    return V
