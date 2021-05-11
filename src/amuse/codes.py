"""AMUSE codes
"""
import sys

_CODES = [
    'athena', 'capreole', 'cachedse',     'gadget2',      'mesa',
    'octgrav',      'twobody', 'capreole',     'hermite0',     'mocassin',
    'phiGRAPE', 'athena',       'evtwin',       'hop',          'seba',
    'bhtree',       'evtwin2sse',   'interface',    'smallN', 'bse',
    'fi',           'mercury',      'sse',
]

__all__ = []


def _import_modules():

    for code in _CODES:
        modulename = 'amuse.legacy.' + code + '.interface'
        try:
            __import__(modulename)
            globals()[code] = sys.modules[modulename]
            __all__.append(code)
        except ImportError as ex:
            modulename = 'amuse.legacy.' + code + '.' + code
            try:
                __import__(modulename) 
                globals()[code] = sys.modules[modulename]
                __all__.append(code)
            except ImportError as ex:
                pass


_import_modules()
