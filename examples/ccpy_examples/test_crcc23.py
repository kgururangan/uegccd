import numpy as np
from hegutils import setup_ueg_driver

def test_crcc23():

    driver = setup_ueg_driver()

    driver.system.print_info()
    driver.run_cc(method="ccsd")
    driver.run_ccp3(method="ccsd(t)")
    driver.run_hbar(method="ccsd")
    driver.run_leftcc(method="left_ccsd")
    driver.run_ccp3(method="crcc23")

    #
    # Check the results
    #
    assert np.allclose(driver.correlation_energy, -0.3061428068, atol=1.0e-08)

if __name__ == "__main__":
    test_crcc23()
