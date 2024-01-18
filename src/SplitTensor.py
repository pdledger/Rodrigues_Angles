import numpy as np

def SplitTensor(TensorArray,Frequencies,N0):
    N=len(Frequencies)
    # Split in to R, I and Rtilde
    Rstore = np.zeros((N,3,3))
    Istore = np.zeros((N,3,3))
    Rtildestore = np.zeros((N,3,3))
    for n in range(N):
        Mlist = TensorArray[n,:]
        Mten = np.array([[Mlist[0], Mlist[1], Mlist[2]],[Mlist[3], Mlist[4], Mlist[5]],[Mlist[6], Mlist[7], Mlist[8]]])
        Rtilde = np.real(Mten)
        I = np.imag(Mten)
        R = np.real(Mten)-N0
        for i in range(3):
            for j in range(3):
                Rstore[n,i,j]=R[i,j]
                Istore[n,i,j]=I[i,j]
                Rtildestore[n,i,j]=Rtilde[i,j]

    return Rstore,Istore,Rtildestore
