# from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import GasolineInterface
# from .interface import Gasoline


class GasolineInterfaceTests(TestWithMPI):

    def test1(self):
        instance = GasolineInterface()
        result, error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
