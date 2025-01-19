import numpy as np
import scipy
import jax
import sympy as sym
import time

def EigenValueDecomposition(N0,TensorArray,Frequencies):

    N = len(Frequencies)
    MultRstore=np.zeros(N)
    MultIstore=np.zeros(N)
    MultRtildestore=np.zeros(N)
    MultN0store=np.zeros(N) # repeated storage of N0 multplicities, eigenvectors, eigenvalues etc

    URstore=np.zeros((N,3),dtype=np.longdouble)
    QRstore=np.zeros((N,3,3),dtype=np.longdouble)
    UIstore=np.zeros((N,3),dtype=np.longdouble)
    QIstore=np.zeros((N,3,3),dtype=np.longdouble)
    UN0store=np.zeros((N,3),dtype=np.longdouble)
    QN0store=np.zeros((N,3,3),dtype=np.longdouble)
    URtildestore=np.zeros((N,3),dtype=np.longdouble)
    QRtildestore=np.zeros((N,3,3),dtype=np.longdouble)
    tol=1e-4
    for n in range(N):
        #print("doing ",n,"of ",N)
        Mlist = TensorArray[n,:]
        Mten = np.array([[Mlist[0], Mlist[1], Mlist[2]],[Mlist[3], Mlist[4], Mlist[5]],[Mlist[6], Mlist[7], Mlist[8]]])#,dtype=np.clongdouble)
        Rtilde = np.real(Mten)#+np.diag(np.random.rand(3))*tol
        I = np.imag(Mten)#+np.diag(np.random.rand(3))*tol
        R = (np.real(Mten)-N0)

        # Computation of eigenvalues, eigenvector
        # Note that the we need to apply to normal matrices and so consider R^TR = R R = R^2
        # Note that the sqrt( eigenvalues of R^2) are the abs(lambda(R))
        # The eigenvectors of R^2 are the same as the eigenvector of R
        #uR,VR = np.linalg.eig((R@R).astype(dtype=float))
        #uRtilde,VRtilde = np.linalg.eig((Rtilde@Rtilde).astype(dtype=float))
        #uI,VI = np.linalg.eig((I@I).astype(dtype=float))
        #uN0,VN0 = np.linalg.eig((N0@N0).astype(dtype=float))
        #uR=np.sqrt(uR)
        #uRtilde=np.sqrt(uRtilde)
        #uI=np.sqrt(uI)
        #uN0=np.sqrt(uN0)
        #print("R")
        uR,VR = np.linalg.eigh(R.astype(dtype=float))#jax.numpy.linalg.eig(R)
        uR=np.real(uR)
        VR=np.real(VR)
        # myR=R.astype(dtype=float)
        # # Create Sympy matrix
        # Rsym=sym.Matrix([[myR[0,0],myR[0,1],myR[0,2]],[myR[1,0],myR[1,1],myR[1,2]],[myR[2,0],myR[2,1],myR[2,2]]]).applyfunc(sym.nsimplify)
        # Rsym=(Rsym+Rsym.T)/2
        # #print(Rsym)
        # out=Rsym.eigenvects()
        # #print(out)
        # for i in range(3):
        #     #(eigenval, multiplicity, eigenspace)
        #     #print(out[i][0],sym.N(sym.re(out[i][0])))
        #     uR[i]=sym.N(sym.re(out[i][0]))
        #     if out[i][1]> 1:
        #         print("multiplicity > 1")
        #     #print(out[i][2][0])
        #     #print((sym.N(out[i][2][0])))
        #     #print(sym.re(sym.N(out[i][2][0])))
        #     #print(np.array(sym.re(sym.N(out[i][2][0]))))
        #     VR[:,i]=(np.array(sym.re(sym.N(out[i][2][0])))).astype(np.float64)[:,0]
        # # make orthonormal
        # VR,dum=np.linalg.qr(VR)


        # for i in range(3):
        #     out=scipy.linalg.null_space(R.astype(dtype=float)-uR[i]*np.eye(3))
        #     n,m=np.shape(out)
        #     if m > 1:
        #         print("Mult more than 1 VR",m)
        #     elif m==1:
        #         VR[:,i]=out[:,0]/np.linalg.norm(out[:,0])
        #     else:
        #         print("Warning VR")


        #print(uR,VR)
        #print("Rtilde")
        uRtilde,VRtilde = np.linalg.eigh(Rtilde.astype(dtype=float))#jax.numpy.linalg.eig(Rtilde)
        uRtilde=np.real(uRtilde)
        VRtilde=np.real(VRtilde)
        #print(uRtilde,Rtilde)
        # for i in range(3):
        #     out=scipy.linalg.null_space(Rtilde.astype(dtype=float)-uRtilde[i]*np.eye(3))
        #     n,m=np.shape(out)
        #     if m > 1:
        #         print("Mult more than 1 VR",m)
        #     elif m==1:
        #         VRtilde[:,i]=out[:,0]/np.linalg.norm(out[:,0])
        #     else:
        #         print("Warning VRtilde")
        # myR=Rtilde.astype(dtype=float)
        # Rtildesym=sym.Matrix([[myR[0,0],myR[0,1],myR[0,2]],[myR[1,0],myR[1,1],myR[1,2]],[myR[2,0],myR[2,1],myR[2,2]]]).applyfunc(sym.nsimplify)
        # Rtildesym=(Rtildesym+Rtildesym.T)/2
        # out=Rtildesym.eigenvects()
        # for i in range(3):
        #     #(eigenval, multiplicity, eigenspace)
        #     #uRtilde[i]=out[i][0]
        #     uRtilde[i]=sym.N(sym.re(out[i][0]))
        #     if out[i][1]> 1:
        #         print("multiplicity > 1")
        #     #VRtilde[:,i]=np.array(out[i][2]).astype(np.float64)[0,:,0]
        #     VRtilde[:,i]=(np.array(sym.re(sym.N(out[i][2][0])))).astype(np.float64)[:,0]
        # # make orthonormal
        # VRtilde,dum=np.linalg.qr(VRtilde)

        #print("I")
        uI,VI = np.linalg.eigh(I.astype(dtype=float))#jax.numpy.linalg.eig(I)
        uI=np.real(uI)
        VI=np.real(VI)
#        for i in range(3):
#            out=scipy.linalg.null_space(I.astype(dtype=float)-uI[i]*np.eye(3))
#            n,m=np.shape(out)
#            if m > 1:
#                print("Mult more than 1 VR",m)
#            elif m==1:
#                print(out[:,0])
#                VI[:,i]=out[:,0]/np.linalg.norm(out[:,0])
#            else:
#                print("Warning VI")

        # myI=I.astype(dtype=float)
        # Isym=sym.Matrix([[myI[0,0],myI[0,1],myI[0,2]],[myI[1,0],myI[1,1],myI[1,2]],[myI[2,0],myI[2,1],myI[2,2]]]).applyfunc(sym.nsimplify)
        # Isym=(Isym+Isym)/2
        # out=Isym.eigenvects()
        # for i in range(3):
        #     #(eigenval, multiplicity, eigenspace)
        #     #uI[i]=out[i][0]
        #     uI[i]=sym.N(sym.re(out[i][0]))
        #     if out[i][1]> 1:
        #         print("multiplicity > 1")
        #     #VI[:,i]=np.array(out[i][2]).astype(np.float64)[0,:,0]
        #     VI[:,i]=(np.array(sym.re(sym.N(out[i][2][0])))).astype(np.float64)[:,0]
        # # make orthonormal
        # VI,dum=np.linalg.qr(VI)




        #print(uI)
        uN0,VN0 = np.linalg.eigh(N0.astype(dtype=float))#jax.numpy.linalg.eig(N0)
        uN0=np.real(uN0)
        VN0=np.real(VN0)
        #uR=np.abs(uR)
        #uRtilde=np.abs(uRtilde)
        #uI=np.abs(uI)
        #uN0=np.abs(uN0)

        # Sort the eigenvalues from largest to smallest magnitude
        # ind = np.argsort(-np.abs(uR))
        # uR = uR[ind]
        # VR = VR[:,ind]
        # ind = np.argsort(-np.abs(uI))
        # uI = uI[ind]
        # VI = VI[:,ind]
        # ind = np.argsort(-np.abs(uRtilde))
        # uRtilde = uRtilde[ind]
        # VRtilde = VRtilde[:,ind]
        # ind = np.argsort(-np.abs(uN0))
        # uN0 = uN0[ind]
        # VN0 = VN0[:,ind]
        # Sort the eigenvalues from largest to smallest

#        ind = np.argsort(-(uR))
#        uR = uR[ind]
#        VR = VR[:,ind]

#        ind = np.argsort(-(uI))
#        uI = uI[ind]
#        VI = VI[:,ind]
#        ind = np.argsort(-(uRtilde))
#        uRtilde = uRtilde[ind]
#        VRtilde = VRtilde[:,ind]
#        ind = np.argsort(-uN0)
#        uN0 = uN0[ind]
#        VN0 = VN0[:,ind]


    # Possible multiplicities are 1, 2 or 3
    # If they are order in terms of increasing multiplicities we will either have
    #  case 1 : 1 , 1, 1 if they are distinct
    #  case 2 : or 1, 2
    #  case 3 : or 3



        MultR = 1#CheckMult(uR,R.astype(dtype=float))
        MultI = 1#CheckMult(uI,I.astype(dtype=float))
        MultRtilde = 1#CheckMult(uRtilde,Rtilde.astype(dtype=float))
        MultN0 = 1#CheckMult(uN0,N0.astype(dtype=float))

        # check ordering of eigenvalues and recompute if two are close
        #uR,VR=checkeigen(R,uR,VR,Frequencies[n])
        #uRtilde,VRtilde=checkeigen(Rtilde,uRtilde,VRtilde,Frequencies[n])
        uI,VI=checkeigen(I,uI,VI,Frequencies[n])
        #print(Frequencies[n],uI,VI)
        #uN0,VN0=checkeigen(N0,uN0,VN0)




        #if n > 5 and n < 8:
        #    print(n)
        #    print("R",uR,VR,MultR)
        #    print("Rtilde",uRtilde,VRtilde,MultRtilde)
        #    print("I",uI,VI,MultI)


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
    Tol=1e-4*np.min(np.abs(u))
    mult=0
#    # Determine the multplicity of eigenvalue lambda_i as 3 - rank(R -lambda_i eye(3))
    for i in range(3):
        mult=np.max([mult,3-np.linalg.matrix_rank(Tensor-u[i]*np.eye(3),tol=Tol)])
#
    # Tol=5e-4
    # mult=1
    # for i in range(3):
    #     for j in range(i+1,3):
    #         if i != j:
    #             if np.abs(u[i]-u[j])/abs(u[i]) < Tol:
    #                 mult=mult+1
    #                 if mult > 3:
    #                     mult =3

    return mult


def checkeigen(M,u,Q,omega):
    order=np.argsort(np.abs(u))
    # Order in absolute ascending order
    u=u[order]
    Q=Q[:,order]
    mult=-1*np.ones(3)
    tag=0
    tol=2e-2#5e-1,1e-1,5e-2,1e-3#5e-2
    for i in range(3):
        for j in range(i+1,3):
            if np.abs(u[i]-u[j])/np.abs(u[i]) < tol and mult[j]==-1:
                mult[i]=j
                mult[j]=i
                tag=1

    if tag==1:
        print("getting new eigenvectors",omega)
        #print(u)
        #print(Q)
        for i in range(3):
            if mult[i]==-1:
                break
        # i contains a first non multuple eigenvalue
        u1=u[i]
        if i==0:
            u2=u[1]
            u3=u[2]
        elif i==1:
            u2=u[0]
            u3=u[2]
        else:
            u2=u[0]
            u3=u[1]

        q1=getevectors(M,u1)
        U,V=orthogcomp(q1)
        q2=computeeigen2(q1,U,V,M,u1,u2)
        q3=np.cross(q1,q2)
        #q2,q3=getevectorsclose(q1)
        Q=np.zeros((3,3))
        Q[:,0]=q1[:]
        Q[:,1]=q2[:]
        Q[:,2]=q3[:]
        u=np.array([u1,u2,u3])
        #print(u)
        #print(Q)
        #time.sleep(10)
    return u,Q
#https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
def getevectors(A,eigenvalue):
    Ameigla=A-eigenvalue*np.eye(3)

    r0=Ameigla[0,:]
    r1=Ameigla[1,:]
    r2=Ameigla[2,:]
    r0xr1=np.cross(r0,r1)
    r0xr2=np.cross(r0,r2)
    r1xr2=np.cross(r1,r2)

    d0=np.dot(r0xr1,r0xr1)
    d1=np.dot(r0xr2,r0xr2)
    d2=np.dot(r1xr2,r1xr2)

    dmax=d0
    imax=0
    if d1 > dmax:
        dmax=d1
        imax=1
    if d2 > dmax:
        imax=2
    if imax==0:
        eigenvector=r0xr1/np.sqrt(d0)
    elif imax==1:
        eigenvector=r0xr2/np.sqrt(d1)
    else:
        eigenvector=r1xr2/np.sqrt(d2)
    return eigenvector


# def getevectorsclose(U):
#     if np.abs(U[0]) > np.abs(U[1]):
#         #The component of max abs value is 0 or 2
#         invlength=1/np.sqrt(U[0]*U[0]+U[2]*U[2])
#         V=np.array([-U[2]*invlength,0,U[0]*invlength])
#     else:
#         #The component of max abs value is 1 or 2
#         invlength=1/np.sqrt(U[1]*U[1]+U[2]*U[2])
#         V=np.array([0,U[2]*invlength,-U[1]*invlength])
#     W=np.cross(U,V)
#     return V,W
def orthogcomp(W):
    if np.abs(W[0]) > np.abs(W[1]):
        #The component of max abs value is 0 or 2
        invlength=1/np.sqrt(W[0]*W[0]+W[2]*W[2])
        U=np.array([-W[2]*invlength,0,W[0]*invlength])
    else:
        #The component of max abs value is 1 or 2
        invlength=1/np.sqrt(W[1]*W[1]+W[2]*W[2])
        U=np.array([0,W[2]*invlength,-W[1]*invlength])
    V=np.cross(W,U)
    return U,V

def computeeigen2(eigen0,U,V,A,eig0,eig1):
    AU = A@U
    AV = A@V
    m00=np.dot(U,AU)-eig1
    m01=np.dot(U,AV)-eig1
    m11=np.dot(V,AV)-eig1
    absM00=np.abs(m00)
    absM01=np.abs(m01)
    absM11=np.abs(m11)
    if absM00>= absM11:
        maxAbsComp=np.max([absM00,absM01])
        if maxAbsComp > 0:
            if absM00 >=absM01:
                m01=m01/m00
                m00=1/np.sqrt(1+m01*m01)
                m01=m01*m00
            else:
                m00=m00/m01
                m01=1/np.sqrt(1+m00*m00)
                m00=m00*m01
            eigenvector1=m01*U-m00*V
        else:
            eigenvector1=U
    else:
        maxAbsComp=np.max([absM11,absM01])
        if maxAbsComp > 0:
            if absM11 >=absM01:
                m01 = m01/m11
                m11=1/np.sqrt(1+m01*m01)
                m01=m01*m11
            else:
                m11=m11/m01
                m01=1/np.sqrt(1+m11*m11)
                m11=m11*m01
            eigenvector1=m11*U - m01*V
        else:
            eigenvector1=U
    return eigenvector1
