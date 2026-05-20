import unittest
from astropy import units as u
import numpy as np
from astropy.io import ascii
from Data_management import column_creator as cc



table_name = "260520_11.44_ST_Catalog_mass_loss_for_0.1-10.0_Gyr_eta-0.1_Rxuv-1.0.ecsv"
catalog = ascii.read(f"Tables/{table_name}")

class TestList(unittest.TestCase):
    def test_add_escape_velocity(self):
        catalog_with_escape_velocity = cc.add_escape_velocity(catalog, "pl_masse", "pl_rade")

        earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
        v_earth = catalog_with_escape_velocity["v_esc"].quantity[earth_idx]
        self.assertAlmostEqual(v_earth, 11.186 * u.km / u.s, places=1)



if __name__ == "__main__":
    unittest.main()