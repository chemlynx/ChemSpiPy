# -*- coding: utf-8 -*-
"""
ChemSpiPy
~~~~~~~~~

Python wrapper for the ChemSpider API.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division


__author__ = 'Matt Swain'
__email__ = 'm.swain@me.com'
__version__ = '1.0.5'
__license__ = 'MIT'


from .api import ChemSpider, MOL2D, MOL3D, BOTH, ASCENDING, DESCENDING, CSID, MASS_DEFECT, MOLECULAR_WEIGHT
from .api import REFERENCE_COUNT, DATASOURCE_COUNT, PUBMED_COUNT, RSC_COUNT
from .objects import Compound
from .search import Results
