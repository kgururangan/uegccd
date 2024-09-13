import numpy as np

def init_heg():
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
