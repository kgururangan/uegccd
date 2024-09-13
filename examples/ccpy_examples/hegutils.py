import numpy as np
from ccpy import Driver
from ccpy.models.integrals import Integral
from ccpy.models.system import System

def init_heg():
    '''Reads basic information about the UEG calcualtion.'''
    with open("ueg.inp", "r") as f:
        for n, line in enumerate(f.readlines()):
            L = line.split()
            if n == 0:
                nelectron, nbasis, nocc = [int(x) for x in L]
            elif n == 1:
                e_hf = float(L[0])
    return nelectron, nbasis, nocc, e_hf

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

def setup_ueg_driver():
    '''Construct a CCpy driver instance populated with the UEG system information and
    bare Hamiltonian.'''
    nelectrons, nbasis, nocc, hf_energy = init_heg()
    nfrozen = 0
    multiplicity = 1

    print("   Reading Fock matrix")
    fock = read_fock(nbasis)

    print("   Reading ERI matrix")
    eri, madelung = read_twobody(nbasis)
    print("")

    system = System(nelectrons,
                    nbasis,
                    multiplicity,
                    nfrozen,
                    reference_energy=hf_energy,
                    nuclear_repulsion=madelung,
                    mo_energies=np.diagonal(fock))

    hamiltonian = Integral(system, 2, {"a": fock, "b": fock,
                                       "aa": eri - eri.transpose(0, 1, 3, 2),
                                       "ab": eri,
                                       "bb": eri - eri.transpose(0, 1, 3, 2)})

    return Driver(system, hamiltonian)

