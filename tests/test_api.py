# -*- coding: utf-8 -*-
"""
test_api
~~~~~~~~

cs = ChemSpider(apikey='woSrjoABv7dpU7C2h3xIUqXHzGTsBWbT')
Test the core API functionality.
response = cs.request('https://api.rsc.org', 'compounds', 'v1', 'lookups/datasources', 'GET')

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
import os
import re

import pytest
import requests
import six
from json import loads

from chemspipy import ChemSpider, MOL2D, MOL3D, BOTH
from chemspipy.errors import ChemSpiPyAuthError, ChemSpiPyServerError


logging.basicConfig(level=logging.WARN)
logging.getLogger('chemspipy').setLevel(logging.DEBUG)

# Security token is retrieved from environment variables
#  CHEMSPIDER_SECURITY_TOKEN = os.environ['CHEMSPIDER_SECURITY_TOKEN']
RSC_APIKEY = os.environ['NEW_CHEMSPIDER_SECURITY_TOKEN']
# Chemspider instances with and without a security token
cs = ChemSpider(apikey=RSC_APIKEY)
cs2 = ChemSpider()


def test_no_security_token():
    """Test ChemSpider can be initialized with no parameters."""
    assert cs2.apikey == None


def test_security_token():
    """Test security token is set correctly when initializing ChemSpider"""
    assert cs.apikey == RSC_APIKEY


def test_chemspider_repr():
    """Test ChemSpider object repr."""
    assert repr(cs) == 'ChemSpider()'
    assert repr(cs2) == 'ChemSpider()'


# MassSpecAPI

def test_get_databases():
    """Test get_databases returns the list of ChemSpider data sources."""
    dbs = cs.get_databases()
    assert all(source in dbs for source in ['Wikipedia', 'ZINC', 'PubChem'])


def test_get_extended_compound_info():
    """Test get_extended_compound_info returns info fo*r a CSID."""
    info = cs.get_extended_compound_info(6543, fields=['smiles', 'commonName','Formula', 'AverageMass', 'MolecularWeight', 'MonoisotopicMass', 'NominalMass', 'ReferenceCount', 'DataSourceCount', 'PubMedCount', 'RSCCount', 'Mol2D', 'Mol3D'])
    assert all(field in info.keys() for field in [
        'id', 'formula', 'smiles', 'averageMass', 'molecularWeight',
        'monoisotopicMass', 'nominalMass', 'commonName'])  #'inchi', 'inchikey',
    assert all(isinstance(info[field], float) for field in [
        'averageMass', 'molecularWeight', 'monoisotopicMass'
    ]) # , 'nominalMass'
    assert isinstance(info['id'], int)
    assert all(isinstance(info[field], six.text_type) for field in [
        'formula', 'smiles', 'commonName'
    ])  #, 'inchi', 'inchikey'

'''
def test_get_extended_compound_info_list():
    """Test get_extended_compound_info_list returns info for a list of CSIDs."""
    info = cs.get_extended_compound_info_list([6543, 1235, 6084])
    assert len(info) == 3
    assert all(field in info[0] for field in [
        'csid', 'molecular_formula', 'smiles', 'inchi', 'inchikey', 'average_mass', 'molecular_weight',
        'monoisotopic_mass', 'nominal_mass', 'alogp', 'xlogp', 'common_name'
    ])
    assert all(isinstance(info[0][field], float) for field in [
        'average_mass', 'molecular_weight', 'monoisotopic_mass', 'nominal_mass', 'alogp', 'xlogp'
    ])
    assert isinstance(info[0]['csid'], int)
    assert all(isinstance(info[0][field], six.text_type) for field in [
        'molecular_formula', 'smiles', 'inchi', 'inchikey', 'common_name'
    ])


def test_get_extended_mol_compound_info_list():
    """Test get_extended_mol_compound_info_list returns info for a list of CSIDs."""
    info = cs.get_extended_mol_compound_info_list([1236], include_external_references=True,
                                                  include_reference_counts=True)
    assert len(info) == 1
    assert all(field in info[0] for field in [
        'csid', 'molecular_formula', 'smiles', 'inchi', 'inchikey', 'average_mass', 'molecular_weight',
        'monoisotopic_mass', 'nominal_mass', 'alogp', 'xlogp', 'common_name', 'reference_count', 'datasource_count',
        'mol_2d'
    ])
    assert all(isinstance(info[0][field], float) for field in [
        'average_mass', 'molecular_weight', 'monoisotopic_mass', 'nominal_mass', 'alogp', 'xlogp'
    ])
    assert all(isinstance(info[0][field], int) for field in ['csid', 'reference_count', 'datasource_count'])
    assert all(isinstance(info[0][field], six.text_type) for field in [
        'molecular_formula', 'smiles', 'inchi', 'inchikey', 'common_name', 'mol_2d'
    ])


def test_get_extended_mol_compound_info_list_dimensions():
    """Test get_extended_mol_compound_info_list returns 2D/3D/both MOL."""
    info = cs.get_extended_mol_compound_info_list([1236], mol_type=MOL2D)
    assert 'mol_2d' in info[0]
    info = cs.get_extended_mol_compound_info_list([1236], mol_type=MOL3D)
    assert 'mol_3d' in info[0]
    info = cs.get_extended_mol_compound_info_list([1236], mol_type=BOTH)
    assert 'mol_2d' in info[0]
    assert 'mol_3d' in info[0]
'''

def test_get_record_mol():
    """Test get_record_mol returns a MOL file."""
    mol = cs.get_record_mol(6084, calc3d=False)
    assert 'V2000' in mol
    assert 'M  END' in mol


def test_simple_search_by_formula():
    """Test simple_search_by_formula returns a list of CSIDs."""
    assert [c.csid for c in cs.simple_search_by_formula('C2H6')] == [6084]


def test_simple_search_by_mass():
    """Test simple_search_by_mass returns a list of CSIDs."""
    csids = [c.csid for c in cs.simple_search_by_mass(17, 0.1)]
    assert len(csids) > 8


# Search

def test_async_simple_search():
    """Test async_simple_search returns a transaction ID."""
    rid = cs.async_simple_search('benzene')
    assert re.compile(r'[a-f0-9\-]{20,50}').search(rid)


def test_async_simple_search_ordered():
    """Test async_simple_search returns a transaction ID."""
    rid = cs.async_simple_search_ordered('glucose')
    assert re.compile(r'[a-f0-9\-]{20,50}').search(rid)


def test_get_async_search_status():
    """Test get_async_search_status returns the status for a transaction ID."""
    rid = cs.async_simple_search('benzene')
    rid = rid['queryId']
    response = cs.get_async_search_status(rid)
    status = response['status']
    assert status in {'Complete', 'Unknown', 'Created', 'Scheduled', 'Processing', 'Suspended', 'PartialResultReady', 'ResultReady'}


def test_get_async_search_status_and_count():
    """Test get_async_search_status_and_count returns the status for a transaction ID."""
    rid = cs.async_simple_search('benzene')
    rid = rid['queryId']
    print(rid)
    print(rid)
    while True:
        status = cs.get_async_search_status_and_count(rid)
        print(status)
        if status['status'] in {'Created', 'Scheduled', 'Processing', 'Complete'}:
            continue
        assert status['count'] == 1
        assert status['message'] == 'Found by approved synonym'
        break


def test_get_async_search_result():
    """Test get_async_search_result returns a list of CSIDs."""
    rid = cs.async_simple_search('benzene')
    while True:
        status = cs.get_async_search_status(rid)
        if status in {'Created', 'Scheduled', 'Processing'}:
            continue
        assert [c.csid for c in cs.get_async_search_result(rid)] == [236]
        break


def test_get_async_search_result_part():
    """Test get_async_search_result_part returns a list of CSIDs."""
    rid = cs.async_simple_search('glucose')
    while True:
        status = cs.get_async_search_status(rid)
        if status in {'Created', 'Scheduled', 'Processing'}:
            continue
        assert len(cs.get_async_search_result_part(rid)) > 6
        assert len(cs.get_async_search_result_part(rid, start=2)) > 2
        assert len(cs.get_async_search_result_part(rid, start=2, count=2)) == 2
        assert len(cs.get_async_search_result_part(rid, start=2, count=99)) > 2
        break


def test_get_compound_info():
    """Test get_compound_info returns info for a CSID."""
    info = cs.get_compound_info(263)
    assert all(field in info for field in ['csid', 'smiles', 'inchi', 'inchikey'])
    assert isinstance(info['csid'], int)
    assert all(isinstance(info[field], six.text_type) for field in ['smiles', 'inchi', 'inchikey'])


def test_get_compound_thumbnail():
    """Test get_compound_thumbnail returns image data for a CSID."""
    img = cs.get_compound_thumbnail(263)
    assert img[:8] == b'\x89PNG\x0d\x0a\x1a\x0a'  # PNG magic number


def test_simple_search():
    """Test simple_search returns a list of CSIDs."""
    assert all(csid in [c.csid for c in cs.simple_search('glucose')] for csid in [5589, 58238, 71358, 96749, 9312824, 9484839])


# Spectra

# This is slow...
# def test_get_all_spectra_info():
#     """Test get_all_spectra_info returns all spectra info."""
#     spectra = cs.get_all_spectra_info()
#     ok_(len(spectra) > 8000)
#     ok_('spectrum_id' in spectrum for spectrum in spectra)

'''
def test_get_spectrum_info():
    """Test get_spectrum_info returns info for the given spectrum ID."""
    info = cs.get_spectrum_info(36)
    assert info['spectrum_id'] == 36
    assert info['csid'] == 235
    assert info['spectrum_type'] == 'HNMR'
    assert info['file_name'] == 'BenzaldehydeHNMR.jdx'
    assert info['submitted_date'] == '2007-08-08T20:18:36.593'


def test_get_compound_spectra_info():
    """Test get_compound_spectra_info returns list of spectra info for the given ChemSpider ID."""
    for s in cs.get_compound_spectra_info(2157):
        assert isinstance(s, dict)
        assert s['csid'] == 2157
        assert isinstance(s['spectrum_id'], int)


def test_get_spectra_info_list():
    """Test get_spectra_info_list returns list of spectra info for a list of CSIDs."""
    assert cs.get_spectra_info_list([263]) == []  # No spectra for this compound
    for s in cs.get_spectra_info_list([2157, 6084]):
        assert s['csid'] in [2157, 6084]
        assert isinstance(s['spectrum_id'], int)
'''

# InChI

def test_get_original_mol():
    """Test get_original_mol returns a MOL file."""
    mol = cs.get_original_mol(6084)
    assert 'V2000' in mol
    assert 'M  END' in mol


# Misc

def test_construct_api_url():
    """Test construction of API URLs."""
    url = cs.construct_api_url('MassSpecAPI', 'GetExtendedCompoundInfo', csid=2157)
    response = requests.get(url)
    assert response.status_code == 200


# Errors

def test_token_needed():
    """Test ChemSpiPyAuthError is raised for certain endpoints if no security_token provided."""
    with pytest.raises(ChemSpiPyAuthError):
        cs2.get_extended_compound_info(263)


def test_invalid_token():
    """Test ChemSpiPyAuthError is raised if a token with invalid format is used."""
    with pytest.raises(ChemSpiPyAuthError):
        mf = ChemSpider('abcde1-1346fa-934a').get_compound(2157).molecular_formula


def test_invalid_token2():
    """Test ChemSpiPyAuthError is raised if a fake token with correct format is used."""
    with pytest.raises(ChemSpiPyAuthError):
        mf = ChemSpider('a1e22457-c835-1234-b141-347bf12fa31c').get_compound(2157).molecular_formula


def test_invalid_rid():
    """Test ChemSpiPyServerError is raised when an invalid transaction ID is used."""
    with pytest.raises(ChemSpiPyServerError):
        cs.get_async_search_status('xxxxxx')


def test_expired_rid():
    """Test ChemSpiPyServerError is raised when a valid but expired transaction ID is used."""
    with pytest.raises(ChemSpiPyServerError):
        cs.get_async_search_status('1a93ee87-acbe-4caa-bc3b-23c3ff39be0f')


def test_fictional_rid():
    """Test ChemSpiPyServerError is raised when a valid but made up transaction ID is used."""
    with pytest.raises(ChemSpiPyServerError):
        cs.get_async_search_status('1a93ee87-acbe-4caa-bc3b-23c3ff39be0a')
