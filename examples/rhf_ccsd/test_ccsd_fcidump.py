import numpy as np
from ccpy import Driver

def test_cct3_hfhminus():
    driver = Driver.from_fcidump(
        fcidump="FCIDUMP",
        nfrozen=0,
        charge=0,
    )
    driver.system.print_info()

    driver.run_cc(method="ccsd")

if __name__ == "__main__":
    test_cct3_hfhminus()
