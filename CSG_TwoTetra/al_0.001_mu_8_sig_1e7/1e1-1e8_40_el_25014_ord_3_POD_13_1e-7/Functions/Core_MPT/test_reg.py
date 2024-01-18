# James Elgy - 04/10/2023

import numpy as np
from matplotlib import pyplot as plt
from ngsolve import *
import scipy.sparse as sp

def test_reg(Theta0Sols, fes2, Theta1, bilinear_bonus_int_order, xivec, inout, epsi, alpha, PODArray, PODTensors, Integration_Order,
             mu_inv, sigma, N0):
    mesh = fes2.mesh

    u, v = fes2.TnT()

    M = BilinearForm(fes2)
    M += SymbolicBFI(inout * u * v, bonus_intorder=bilinear_bonus_int_order)
    M.Assemble()

    M_out = BilinearForm(fes2)
    M_out += SymbolicBFI((1-inout) * u * v, bonus_intorder=bilinear_bonus_int_order)
    M_out.Assemble()

    rows, cols, vals = M.mat.COO()
    del M
    M = sp.csr_matrix((vals, (rows, cols)))
    del rows, cols, vals

    rows, cols, vals = M_out.mat.COO()
    del M_out
    M_out = sp.csr_matrix((vals, (rows, cols)))
    del rows, cols, vals

    s1 = LinearForm(fes2)
    s1 += SymbolicLFI(inout * xivec[0] * v)
    s2 = LinearForm(fes2)
    s2 += SymbolicLFI(inout * xivec[1] * v)
    s3 = LinearForm(fes2)
    s3 += SymbolicLFI(inout * xivec[2] * v)

    s1_out = LinearForm(fes2)
    s1_out += SymbolicLFI((1-inout) * xivec[0] * v)
    s2_out = LinearForm(fes2)
    s2_out += SymbolicLFI((1-inout) * xivec[1] * v)
    s3_out = LinearForm(fes2)
    s3_out += SymbolicLFI((1-inout) * xivec[2] * v)

    s1.Assemble()
    s2.Assemble()
    s3.Assemble()
    s1_out.Assemble()
    s2_out.Assemble()
    s3_out.Assemble()


    term1 = np.zeros((Theta1.shape[1], 3, 3), dtype=complex)
    term2 = np.zeros((Theta1.shape[1], 3, 3), dtype=complex)
    term3 = np.zeros((Theta1.shape[1], 3, 3), dtype=complex)
    term4 = np.zeros((Theta1.shape[1], 3, 3), dtype=complex)
    for k in range(Theta1.shape[1]):
       
        for i in range(3):
            match i:
                case 0:
                    s = s1.vec.FV().NumPy()[:]
                    s_out = s1_out.vec.FV().NumPy()[:]
                case 1:
                    s = s2.vec.FV().NumPy()[:]
                    s_out = s2_out.vec.FV().NumPy()[:]
                case 2:
                    s = s3.vec.FV().NumPy()[:]
                    s_out = s3_out.vec.FV().NumPy()[:]

            t0 = Theta0Sols[:,i]
            for j in range(3):
                t1j = np.squeeze(Theta1[:, k, j])
                t1i = np.squeeze(Theta1[:, k, i])

                # Term 1:
                part1 = t1j[None,:] @ M @ t0[:,None] # Using [None,:] as quick transpose 
                part2 = s[None,:] @ t1j
                term1[k, i, j] = epsi * (alpha**3)/4 * (part1 + part2)

                # Term 2:
                term2[k, i, j] = epsi * (alpha**3)/4 * np.conj((t1i[None,:]) @ M_out @ t1j[:, None])

                # Term 3:
                term3[k, i, j] = epsi * (alpha**3)/4 * s_out[None,:] @ t1j[:,None]

    # Computing term4 as imag part of R and I:
    Mu0 = 4* np.pi * 10 ** (-7)
    nu_no_omega = Mu0 * (alpha ** 2)
    Theta1i = GridFunction(fes2)
    Theta0i = GridFunction(fes2)
    Theta0j = GridFunction(fes2)
    Theta1j = GridFunction(fes2)
    R_im = np.zeros((3,3), dtype=complex)
    I_im = np.zeros((3,3), dtype=complex)
    for k, omega in enumerate(PODArray):
        for i in range(3):
            Theta1i.vec.FV().NumPy()[:] = Theta1[:, k, i]
            Theta0i.vec.FV().NumPy()[:] = Theta0Sols[:, i]
            xii = xivec[i]
            for j in range(3):
                Theta0j.vec.FV().NumPy()[:] = Theta0Sols[:, j]
                xij = xivec[j]
                Theta1j.vec.FV().NumPy()[:] = Theta1[:, k, j]
                # Real and Imaginary parts
                with TaskManager():
                    R_im[i, j] = -(((alpha ** 3) / 4) * Integrate((mu_inv) * (curl(Theta1j) * Conj(curl(Theta1i))),
                                                                mesh, order=Integration_Order)).imag
                    I_im[i, j] = ((alpha ** 3) / 4) * Integrate(inout * nu_no_omega * omega * sigma * (
                                (Theta1j + Theta0j + xij) * (Conj(Theta1i) + Theta0i + xii)), mesh,
                                                                order=Integration_Order).imag

                    term4[k, i,j] = R_im[i,j] + 1j * I_im[i,j]
    
    term1_norm_real = []
    term2_norm_real = []
    term3_norm_real = []
    term4_norm_real = []
    term1_norm_imag = []
    term2_norm_imag = []
    term3_norm_imag = []
    term4_norm_imag = []
    for k in range(Theta1.shape[1]):
        term1_norm_real += [np.linalg.norm(np.squeeze(term1[k,:,:].real))]
        term2_norm_real += [np.linalg.norm(np.squeeze(term2[k,:,:].real))]
        term3_norm_real += [np.linalg.norm(np.squeeze(term3[k,:,:].real))]
        term4_norm_real += [np.linalg.norm(np.squeeze(term4[k,:,:].real))]

        term1_norm_imag += [np.linalg.norm(np.squeeze(term1[k,:,:].imag))]
        term2_norm_imag += [np.linalg.norm(np.squeeze(term2[k,:,:].imag))]
        term3_norm_imag += [np.linalg.norm(np.squeeze(term3[k,:,:].imag))]
        term4_norm_imag += [np.linalg.norm(np.squeeze(term4[k,:,:].imag))]


    plt.figure()
    plt.semilogx(PODArray, term1_norm_real, label='Term 1', marker='x')
    plt.semilogx(PODArray, term2_norm_real, label='Term 2', marker='^')
    plt.semilogx(PODArray, term3_norm_real, label='Term 3', marker='+')
    plt.semilogx(PODArray, term4_norm_real, label='Term 4', marker='*')
    
    plt.legend()
    plt.xlabel('$\omega$, [rad/s]')
    plt.ylabel('$||Re(y)||_F$')
    plt.yscale('log')

    plt.figure()
    plt.semilogx(PODArray, term1_norm_imag, label='Term 1', marker='x')
    plt.semilogx(PODArray, term2_norm_imag, label='Term 2', marker='^')
    plt.semilogx(PODArray, term3_norm_imag, label='Term 3', marker='+')
    plt.semilogx(PODArray, term4_norm_imag, label='Term 4', marker='*')

    plt.legend()
    plt.yscale('log')
    plt.xlabel('$\omega$, [rad/s]')
    plt.ylabel('$||Im(y)||_F$')
    #plt.show()

    R_norm = []
    R_tilde_norm = []
    I_norm = []
    Z_upper_bound = []
    # norms of R and I
    for k in range(len(PODArray)):
        R = PODTensors[k, :].reshape(3,3).real - N0
        R_tilde = PODTensors[k,:].reshape(3,3).real
        I = PODTensors[k, :].reshape(3,3).imag
        R_norm += [np.linalg.norm(R)]
        I_norm += [np.linalg.norm(I)]
        R_tilde_norm += [np.linalg.norm(R_tilde)]

    R_norm = np.asarray(R_norm)
    R_tilde_norm = np.asarray(R_tilde_norm)
    I_norm = np.asarray(I_norm)
    R_epsi_norm = np.asarray(term1_norm_real) + np.asarray(term2_norm_imag) + np.asarray(term3_norm_real) + np.asarray(term4_norm_real)
    I_epsi_norm = np.asarray(term1_norm_imag) + np.asarray(term2_norm_imag) + np.asarray(term3_norm_imag) + np.asarray(term4_norm_imag)

    for k in range(len(PODArray)):
        Z_upper_bound += [2 *(R_epsi_norm[k] * I_norm[k]) + 2*(I_epsi_norm[k] * R_norm[k]) + 2*(R_epsi_norm[k] * I_epsi_norm[k])]


    plt.figure()
    plt.loglog(PODArray, Z_upper_bound)
    plt.ylabel('$||Z||_F$ upper bound')
    plt.xlabel('$\omega$, [rad/s]')


    
    ## Computing regularisation term for N0:
    Theta0i = GridFunction(fes2)
    Theta0j = GridFunction(fes2)
    N0_epsi = np.zeros((3,3), dtype=complex)
    for i in range(3):
        Theta0i.vec.FV().NumPy()[:] = Theta0Sols[:, i]
        for j in range(3):
            Theta0j.vec.FV().NumPy()[:] = Theta0Sols[:, j]
            with TaskManager():
                N0_epsi[i, j] = -((alpha ** 3) / 4) * Integrate( epsi * Theta0i * Theta0j, mesh, order=Integration_Order)

    
    # Computing upper bounds on Z_tilde
    R_tilde_epsi = np.zeros((Theta1.shape[1], 3,3), dtype=complex)
    R_tilde_epsi_norm = np.zeros((Theta1.shape[1]))
    for k in range(len(PODArray)):
        R_tilde_epsi[k,:,:] = np.asarray(term1[k,:,:]) + np.asarray(term2[k,:,:]) + np.asarray(term3[k,:,:]) + np.asarray(term4[k,:,:]) + N0_epsi
        R_tilde_epsi_norm[k] = np.linalg.norm(R_tilde_epsi[k,:,:])

    Z_tilde_upper_bound = []
    for k in range(len(PODArray)):
        Z_tilde_upper_bound += [2 *(R_tilde_epsi_norm[k] * I_norm[k]) + 2*(I_epsi_norm[k] * R_tilde_norm[k]) + 2*(R_tilde_epsi_norm[k] * I_epsi_norm[k])]

    
    # Checking |R| compared to regularisation terms.
    
    sum_reg_terms_real = np.zeros((len(PODArray), 3, 3), dtype=complex)
    sum_reg_terms_imag = np.zeros((len(PODArray), 3, 3), dtype=complex)
    for k in range(len(PODArray)):
        sum_reg_terms_real[k, :, :] = np.abs(term1[k,:,:].real) + np.abs(term2[k,:,:].real) + np.abs(term3[k,:,:].real)  + np.abs(term4[k,:,:].real) 
        sum_reg_terms_imag[k, :, :] = np.abs(term1[k,:,:].imag) + np.abs(term2[k,:,:].imag) + np.abs(term3[k,:,:].imag)  + np.abs(term4[k,:,:].imag) 

    sum_reg_terms_real = sum_reg_terms_real.reshape(len(PODArray), 9)
    sum_reg_terms_imag = sum_reg_terms_imag.reshape(len(PODArray), 9)
    plt.figure()
    label_list = [12, 13, 21, 23, 31, 32]
    cols = ['tab:blue', 'tab:green', 'tab:red', 'tab:pink', 'tab:orange', 'tab:brown']
    for ind, coeff in enumerate([1,2,3,5,6,7]):
        plt.loglog(PODArray, np.abs(PODTensors[:, coeff].real), label=f'$|\mathcal{{R}}_{{{label_list[ind]}}}|$', color=cols[ind])
        plt.loglog(PODArray, np.abs(sum_reg_terms_real[:, coeff].real), label=f'$|\Delta(\mathcal{{R}})_{{{label_list[ind]}}}|$', linestyle='--', color=cols[ind])
    plt.legend()
    plt.ylabel('Real off diag')
    plt.xlabel('$\omega$, [rad/s]')
    
    plt.figure()
    for ind, coeff in enumerate([1,2,3,5,6,7]):
        plt.loglog(PODArray, np.abs(PODTensors[:, coeff].imag), label=f'$|\mathcal{{I}}_{{{label_list[ind]}}}|$', color=cols[ind])
        plt.loglog(PODArray, np.abs(sum_reg_terms_imag[:, coeff].real), label=f'$|\Delta(\mathcal{{I}})_{{{label_list[ind]}}}|$', linestyle='--', color=cols[ind])
    plt.legend()
    plt.ylabel('Imag off diag')
    plt.xlabel('$\omega$, [rad/s]')
    
    np.savetxt('sum_reg_terms_real.csv', sum_reg_terms_real)
    np.savetxt('sum_reg_terms_imag.csv', sum_reg_terms_imag)
    
    return term1, term2, term3, term4, Z_upper_bound, Z_tilde_upper_bound

if __name__ == '__main__':
    pass
