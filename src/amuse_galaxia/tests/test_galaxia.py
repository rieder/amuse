from amuse.support.testing import amusetest
from amuse.units import units, nbody_system, constants

from amuse_galaxia.interface import Galaxia


class TestGalaxia(amusetest.TestCase):
    def test1(self):
        galaxy = Galaxia()
        self.assertEqual(galaxy.model_time, 0. | units.yr)

    def test2(self):
        galaxy = Galaxia()
        vc1 = galaxy.get_velcirc(10. | units.kpc, 0. | units.kpc, 0. | units.kpc)
        vc2 = galaxy.get_velcirc(0. | units.kpc, 10. | units.kpc, 0. | units.kpc)
        self.assertEqual(vc1, vc2)
        self.assertAlmostEqual(vc1, 219.324003066 | units.kms)
