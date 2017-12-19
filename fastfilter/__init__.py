from pkg_resources import get_distribution

try:
    __version__ = get_distribution('fastfilter').version
except:
    __version__ = 'local'


__all__ = [
    'makedb',
    'runfilter'
]

from fastfilter import *