from amuse.test.amusetest import TestWithMPI

from .interface import MetisseInterface
from .interface import Metisse

class MetisseInterfaceTests(TestWithMPI):
    def test_initialize(self):
        instance = MetisseInterface()
        error = instance.initialize()
        self.assertEqual(error, 0)
        instance.stop()


class MetisseTests(TestWithMPI):
    def test_initialize(self):
        instance = Metisse()
        error = instance.initialize()
        self.assertEqual(error, 0)
        instance.stop()
