import numpy as np
from hegutils import init_heg, read_fock, read_twobody
from miniccpy.printing import print_custom_system_information_rhf
from miniccpy.driver import run_cc_calc

def test_ccsd():

    nelectrons, nbasis, nocc, hf_energy = init_heg()
    nfrozen = 0

    o = slice(0, nocc)
    v = slice(nocc, nbasis)

    print("   Reading Fock matrix")
    fock = read_fock(nbasis)

    print("   Reading ERI matrix")
    g, _ = read_twobody(nbasis)
    print("")

    print_custom_system_information_rhf(fock, nelectrons, nfrozen, hf_energy)

    T, E_corr = run_cc_calc(fock, g, o, v, method='rccsd', maxit=80)

    #
    # Check the results
    #
    assert np.allclose(E_corr, -0.1369034549, atol=1.0e-08) 

if __name__ == "__main__":
    test_ccsd()
