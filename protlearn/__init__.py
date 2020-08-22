from .preprocessing import encode

from .features.ngram import ngram
from .features.posrich import posrich
from .features.entropy import entropy
from .features.paac import paac
from .features.cksaap import cksaap
from .features.ctd import ctd
from .features.binary import binary
from .features.atc import atc
from .features.moran import moran
from .features.geary import geary
from .features.moreau_broto import moreau_broto

from .utils.validation import check_input

__version__ = '1.9.2'

__all__ = ['preprocessing',
           'features',
           'dimreduction']