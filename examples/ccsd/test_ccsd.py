import numpy as np
from miniccpy.printing import print_custom_system_information, print_custom_system_information_rhf
from miniccpy.driver import run_cc_calc, run_mpn_calc
from miniccpy.integrals import spatial_to_spinorb

def init_heg():
    with open("nocc", "r") as f:
        nocc = int(f.readlines()[0])
    with open("nelectrons", "r") as f:
        nelectron = int(f.readlines()[0])
    with open("nbasis", "r") as f:
        norb = int(f.readlines()[0])
    return nelectron, norb, nocc

def read_fock(nbasis):
    '''Reads the Fock matrix in onebody.inp.'''
    fock = np.zeros((nbasis, nbasis))
    with open("onebody.inp", "r") as f:
        lines = f.readlines()
        cnt = 0
        for p in range(nbasis):
            for q in range(p + 1):
                fock[p, q] = float(lines[cnt].split()[0])
                fock[q, p] = fock[p, q]
                cnt += 1
    return fock 

def read_twobody(nbasis):
    '''Reads the ERI matrix (in physics notation) in twobody.inp.'''
    eri = np.zeros((nbasis, nbasis, nbasis, nbasis))
    with open("twobody.inp", "r") as f:
        for line in f:
            L = line.split()
            p, q, r, s = [int(x) - 1 for x in L[:-1]]
            val = float(L[-1])
            # nuclear repulsion value
            if p == -1 and q == -1 and r == -1 and s == -1:
                nuclear_repulsion = val
            else:
                eri[p, q, r, s] = val
    return eri, nuclear_repulsion

def test_ccsd():

    nelectrons, nbasis, _ = init_heg()
    nfrozen = 0
    hf_energy = 0.0

    o = slice(nelectrons)
    v = slice(nelectrons, nbasis * 2)

    print("Reading Fock matrix")
    fock = read_fock(nbasis)
    print("Reading ERI matrix")
    g, _ = read_twobody(nbasis)
    print("Converting to spinorbitals")
    fock, g = spatial_to_spinorb(fock, g)
    g -= np.transpose(g, (0, 1, 3, 2))

    print_custom_system_information(fock, nelectrons, nfrozen, hf_energy)

    e_mp2 = run_mpn_calc(fock, g, o, v, method='mp2')
    T, E_corr = run_cc_calc(fock, g, o, v, method='ccsd', maxit=80)

    #
    # Check the results
    #
    assert np.allclose(E_corr, -0.0907402833, atol=1.0e-08)

if __name__ == "__main__":
    test_ccsd()



