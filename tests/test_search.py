# -*- coding: utf-8 -*-
"""
test_search
~~~~~~~~~~~

Test the search wrapper.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
import os

import pytest

from chemspipy import ChemSpider, ASCENDING, DESCENDING, CSID, REFERENCE_COUNT, MOLECULAR_WEIGHT
from chemspipy.errors import ChemSpiPyServerError


logging.basicConfig(level=logging.WARN, format='%(levelname)s:%(name)s:(%(threadName)-10s):%(message)s')
logging.getLogger('chemspipy').setLevel(logging.DEBUG)

# API key is retrieved from environment variables
CHEMSPIDER_API_KEY = os.environ['CHEMSPIDER_API_KEY']
cs = ChemSpider(CHEMSPIDER_API_KEY)


def test_search_smiles():
    """Test SMILES input to search."""
    results = cs.search('O=C(OCC)C')
    assert results.ready() == False
    results.wait()
    assert results.ready() == True
    assert results.success() == True
    assert results.message == 'Found by conversion query string to chemical structure (full match)'
    assert results[0].csid == 8525
    assert results.duration.total_seconds() > 0


def test_search_csid():
    """Test ChemSpider ID input to search."""
    results = cs.search(8525)
    assert results.message == 'Found by CSID'
    assert len(results) == 1
    assert repr(results) == 'Results([Compound(8525)])'
    assert results[0].csid == 8525


def test_search_name():
    """Test name input to search."""
    results = cs.search('propanol')
    assert results.message == 'Found by approved synonym'
    assert results[0].csid == 1004


def test_search_iter():
    """Test iteration of search results."""
    for result in cs.search('glucose'):
        assert isinstance(result.csid, int)


def test_search_ordered_csid():
    """Test search results ordered by CSID."""
    results = cs.search('glucose', order=CSID)
    assert list(results)  == sorted(results, key=lambda x: x.csid)


def test_search_ordered_csid_descending():
    """Test search results ordered by CSID and direction descending."""
    results = cs.search('glucose', order=CSID, direction=DESCENDING)
    assert list(results) == sorted(results, key=lambda x: x.csid, reverse=True)


def test_search_ordered_ref_descending():
    """Test search results ordered by CSID and direction descending."""
    results = cs.search('glucose', order=REFERENCE_COUNT, direction=DESCENDING)
    assert [result.csid for result in results]


def test_search_ordered_weight_ascending():
    """Test search results ordered by CSID and direction descending."""
    results = cs.search('P', order=MOLECULAR_WEIGHT, direction=ASCENDING)
    assert list(results) == sorted(results, key=lambda x: x.molecular_weight)


def test_search_no_results():
    """Test name input to search."""
    results = cs.search('aergherguyaelrgiaubrfawyef')
    assert results.message == 'No results found'
    assert results.ready() == True
    assert results.success() == True
    assert len(results) == 0


def test_too_high_index():
    """Test IndexError is raised for a too high index."""
    with pytest.raises(IndexError):
        result = cs.search('glucose')[7843]


def test_search_failed():
    """Test ChemSpiPyServerError is raised for an invalid SMILES."""
    results = cs.search('O=C(OCC)C*')
    results.wait()
    assert isinstance(results.exception, ChemSpiPyServerError)
    assert results.status == 'Failed'
    assert repr(results) == 'Results(Failed)'
    assert results.ready() == True
    assert results.success() == False
    assert results.count == 0
    assert results.duration.total_seconds() > 0


def test_search_exception():
    """Test ChemSpiPyServerError is raised for an invalid SMILES."""
    with pytest.raises(ChemSpiPyServerError):
        results = cs.search('O=C(OCC)C*', raise_errors=True)
        results.wait()


# ordered search - ascending/descending, different sort orders
