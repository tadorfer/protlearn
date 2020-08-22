from preprocessing import encode

from features import ngram
from features import posrich
from features import entropy
from features import paac
from features import cksaap
from features import ctd
from features import binary
from features import atc
from features import moran
from features import geary
from features import moreau_broto

from utils.validation import check_input

__version__ = '1.9.2'

__all__ = ['preprocessing',
           'features',
           'dimreduction']