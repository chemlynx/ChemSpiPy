# -*- coding: utf-8 -*-
"""
test_compound
~~~~~~~~~~~~~

Test the Compound object.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
import os

import requests

from chemspipy import ChemSpider, Compound #  , Spectrum


logging.basicConfig(level=logging.WARN, format='%(levelname)s:%(name)s:(%(threadName)-10s):%(message)s')
logging.getLogger('chemspipy').setLevel(logging.DEBUG)

# Security token is retrieved from environment variables
RSC_APIKEY = os.environ['NEW_CHEMSPIDER_SECURITY_TOKEN']
cs = ChemSpider(apikey=RSC_APIKEY)


def test_get_compound():
    """Test getting a compound by ChemSpider ID."""
    compound = cs.get_compound(2157)
    assert isinstance(compound, Compound)
    assert compound.csid == 2157
    compound = cs.get_compound('2157')
    assert isinstance(compound, Compound)
    assert compound.csid == 2157


def test_get_compounds():
    """Test getting multiple compounds by ChemSpider ID."""
    compounds = cs.get_compounds([2157, 13837760])
    assert [c.csid for c in compounds], [2157 == 13837760]
    for c in compounds:
        assert 'http://' in c.image_url
        assert c.average_mass > 0


def test_compound_init():
    """Test instantiating a Compound directly."""
    compound = Compound(cs, 2157)
    assert compound.csid == 2157


def test_compound_equality():
    """Test equality test by ChemSpider ID."""
    c1 = cs.get_compound(13837760)
    c2 = cs.get_compound(2157)
    c3 = cs.get_compound(2157)
    assert c1 != c2
    assert c2 == c3


def test_compound_repr():
    """Test Compound object repr."""
    assert repr(cs.get_compound(1234)) == 'Compound(1234)'


def test_image_url():
    """Test image_url returns a valid URL."""
    url = cs.get_compound(2157).image_url
    response = requests.get(url)
    assert 'http://www.chemspider.com/ImagesHandler.ashx?id=' in url
    assert response.status_code == 200


def test_molecular_formula():
    """Test Compound property molecular_formula."""
    compound = cs.get_compound(2157)
    assert compound.molecular_formula == 'C_{9}H_{8}O_{4}'
    # Ensure value is the same on subsequent access from cache
    assert compound.molecular_formula == 'C_{9}H_{8}O_{4}'


def test_smiles():
    """Test Compound property smiles."""
    compound = cs.get_compound(2157)
    assert compound.smiles == 'CC(=O)Oc1ccccc1C(=O)O'
    # Ensure value is the same on subsequent access from cache
    assert compound.smiles == 'CC(=O)Oc1ccccc1C(=O)O'


def test_inchi():
    """Test Compound property inchi."""
    compound = cs.get_compound(2157)
    assert compound.inchi == 'InChI=1/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
    # Ensure value is the same on subsequent access from cache
    assert compound.inchi == 'InChI=1/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'


def test_stdinchi():
    """Test Compound property stdinchi."""
    compound = cs.get_compound(2157)
    assert compound.stdinchi == 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
    # Ensure value is the same on subsequent access from cache
    assert compound.stdinchi == 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'


def test_inchikey():
    """Test Compound property inchikey."""
    compound = cs.get_compound(2157)
    assert compound.inchikey == 'BSYNRYMUTXBXSQ-UHFFFAOYAW'
    # Ensure value is the same on subsequent access from cache
    assert compound.inchikey == 'BSYNRYMUTXBXSQ-UHFFFAOYAW'


def test_stdinchikey():
    """Test Compound property stdinchikey."""
    compound = cs.get_compound(2157)
    assert compound.stdinchikey == 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
    # Ensure value is the same on subsequent access from cache
    assert compound.stdinchikey == 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'


def test_masses():
    """Test Compound property average_mass, molecular_weight, monoisotopic_mass, nominal_mass."""
    compound = cs.get_compound(2157)
    assert 180 < compound.average_mass < 180.2
    assert 180 < compound.molecular_weight < 180.2
    assert 180 < compound.monoisotopic_mass < 180.2
    assert compound.nominal_mass == 180


def test_descriptors():
    """Test Compound property alogp, xlogp."""
    compound = cs.get_compound(348191)
    assert compound.alogp == 0.0
    assert compound.xlogp == 1.2


def test_name():
    """Test Compound property common_name."""
    compound = cs.get_compound(2157)
    assert compound.common_name == 'Aspirin'


def test_molfiles():
    """Test Compound property mol2d, mol3d, mol_raw."""
    compound = cs.get_compound(2157)
    assert 'V2000' in compound.mol_2d
    assert 'V2000' in compound.mol_3d
    assert 'V2000' in compound.mol_raw


def test_image():
    """Test Compound property image."""
    compound = cs.get_compound(2157)
    assert compound.image[:8] == b'\x89PNG\x0d\x0a\x1a\x0a'  # PNG magic number

'''
def test_spectra():
    """Test Compound property spectra."""
    compound = cs.get_compound(2157)
    for s in compound.spectra:
        assert isinstance(s, Spectrum)
        assert s.csid == 2157
        assert isinstance(s.spectrum_id, int)
    compound = cs.get_compound(263)
    assert compound.spectra == []
'''
