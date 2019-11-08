from amuse.test.amusetest import TestWithMPI

from .interface import PentacleInterface
from .interface import Pentacle

class PentacleInterfaceTests(TestWithMPI):
    
    def xtest1(self):
        instance = PentacleInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()

class PentacleTests(TestWithMPI):

    def test_startstop(self):
        instance = Pentacle()
        parameters = instance.parameters
        print(parameters)
        instance.stop()
