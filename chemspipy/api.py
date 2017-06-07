# -*- coding: utf-8 -*-
"""
chemspipy.api
~~~~~~~~~~~~~

Core API for interacting with ChemSpider web services.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from base64 import b64decode
import logging
import sys
import warnings

try:
    from lxml import etree
except ImportError:
    try:
        import xml.etree.cElementTree as etree
    except ImportError:
        import xml.etree.ElementTree as etree

import requests
import six

from . import __version__
from .errors import ChemSpiPyError, ChemSpiPyParseError, ChemSpiPyAuthError, ChemSpiPyServerError
from .errors import ChemSpiPyNotFoundError
from .objects import Compound, Spectrum
from .search import Results


log = logging.getLogger(__name__)

#: 2D coordinate dimensions
MOL2D = '2d'
#: 3D coordinate dimensions
MOL3D = '3d'
#: Both coordinate dimensions
BOTH = 'both'

#: Ascending sort direction
ASCENDING = 'ascending'
#: Descending sort direction
DESCENDING = 'descending'

#: CSID sort order
CSID = 'csid'
#: Mass defect sort order
MASS_DEFECT = 'mass_defect'
#: Molecular weight sort order
MOLECULAR_WEIGHT = 'molecular_weight'
#: Reference count sort order
REFERENCE_COUNT = 'reference_count'
#: Datasource count sort order
DATASOURCE_COUNT = 'datasource_count'
#: Pubmed count sort order
PUBMED_COUNT = 'pubmed_count'
#: RSC count sort order
RSC_COUNT = 'rsc_count'


#: Coordinate dimensions
DIMENSIONS = {
    MOL2D: 'e2D',
    MOL3D: 'e3D',
    BOTH: 'eBoth'
}

#: Sort directions
DIRECTIONS = {
    ASCENDING: 'eAscending',
    DESCENDING: 'eDescending'
}

#: Sort orders
ORDERS = {
    CSID: 'eCSID',
    MASS_DEFECT: 'eMassDefect',
    MOLECULAR_WEIGHT: 'eMolecularWeight',
    REFERENCE_COUNT: 'eReferenceCount',
    DATASOURCE_COUNT: 'eDataSourceCount',
    PUBMED_COUNT: 'ePubMedCount',
    RSC_COUNT: 'eRscCount'
}

#: API to python field mappings
FIELDS = {
    'CSID': ('csid', int),
    'csid': ('csid', int),
    'MF': ('molecular_formula', six.text_type),
    'SMILES': ('smiles', six.text_type),
    'InChI': ('inchi', six.text_type),
    'InChIKey': ('inchikey', six.text_type),
    'AverageMass': ('average_mass', float),
    'MolecularWeight': ('molecular_weight', float),
    'MonoisotopicMass': ('monoisotopic_mass', float),
    'NominalMass': ('nominal_mass', float),
    'ALogP': ('alogp', float),
    'XLogP': ('xlogp', float),
    'CommonName': ('common_name', six.text_type),
    'MOL2d': ('mol_2d', six.text_type),
    'MOL3d': ('mol_3d', six.text_type),
    'ReferenceCount': ('reference_count', int),
    'DataSourceCount': ('datasource_count', int),
    'PubMedCount': ('pubmed_count', int),
    'RSCCount': ('rsc_count', int),
    'ExternalReferences': ('external_references', list),
    'ds_name': ('datasource_name', six.text_type),
    'ds_url': ('datasource_url', six.text_type),
    'ext_id': ('external_id', six.text_type),
    'ext_url': ('external_url', six.text_type),
    'Status': ('status', six.text_type),
    'Count': ('count', int),
    'Message': ('message', six.text_type),
    'Elapsed': ('elapsed', six.text_type),
    'spc_id': ('spectrum_id', int),
    'spc_type': ('spectrum_type', six.text_type),
    'file_name': ('file_name', six.text_type),
    'comments': ('comments', six.text_type),
    'original_url': ('original_url', six.text_type),
    'submitted_date': ('submitted_date', six.text_type),
}


class BaseChemSpider(object):

    def __init__(self, security_token=None, user_agent=None, api_url=None):
        """

        :param string security_token: (Optional) Your ChemSpider security token.
        :param string user_agent: (Optional) Identify your application to ChemSpider servers.
        :param string api_url: (Optional) Alternative API server.
        """
        log.debug('Initializing ChemSpider')
        self.api_url = api_url if api_url else 'https://www.chemspider.com'
        self.http = requests.session()
        self.http.headers['User-Agent'] = user_agent if user_agent else 'ChemSpiPy/%s Python/%s ' % (__version__, sys.version.split()[0])
        self.security_token = security_token

    def request(self, api, endpoint, **params):
        """Construct API request and return the XML response.

        :param string api: The specific ChemSpider API to call (MassSpec, Search, Spectra, InChI).
        :param string endpoint: ChemSpider API endpoint.
        :param params: (Optional) Parameters for the ChemSpider endpoint as keyword arguments.
        :rtype: xml tree
        """
        url = '%s/%s.asmx/%s' % (self.api_url, api, endpoint)
        log.debug('Request: %s %s', url, params)
        params['token'] = self.security_token
        try:
            response = self.http.post(url, data=params)
        except requests.RequestException as e:
            raise ChemSpiPyError(six.text_type(e))
        if response.status_code == 500:
            if 'Missing parameter: token.' in response.text:
                raise ChemSpiPyAuthError('Endpoint requires a security token.')
            elif 'Error converting data type nvarchar to uniqueidentifier' in response.text:
                # Generally when supplying a security token with incorrect format
                raise ChemSpiPyAuthError('Invalid security token. Did you copy the entire token?')
            elif 'Unauthorized web service usage' in response.text:
                # Fake/incorrect token (but in correct format)
                raise ChemSpiPyAuthError(response.text)
            elif 'Unable to get record details' in response.text:
                # Generally when requesting a non-existent CSID
                raise ChemSpiPyNotFoundError(response.text)
            elif 'Unable to get records spectra' in response.text:
                # No spectra for a CSID, shouldn't be an exception
                return []
            else:
                raise ChemSpiPyServerError(response.text)
        try:
            tree = etree.fromstring(response.content)
        except etree.ParseError as e:
            raise ChemSpiPyParseError('Unable to parse XML response: %s' % e)
        return tree

    def construct_api_url(self, api, endpoint, **params):
        """Construct a Chemspider API url, encoded, with parameters as a GET querystring.

        :param string api: The specific ChemSpider API to call (MassSpecAPI, Search, Spectra, InChI).
        :param string endpoint: ChemSpider API endpoint.
        :param params: (Optional) Parameters for the ChemSpider endpoint as keyword arguments.
        :rtype: string
        """
        querystring = []
        for k, v in params.items():
            querystring.append('%s=%s' % (k, six.moves.urllib.parse.quote_plus(six.text_type(v))))
        if self.security_token:
            querystring.append('token=%s' % self.security_token)
        return '%s/%s.asmx/%s?%s' % (self.api_url, api, endpoint, '&'.join(querystring))


def xml_to_dict(t):
    """Convert a ChemSpider XML response to a python dict."""
    d = {}
    for child in t:
        tag = child.tag.split('}')[1]
        tag, rtype = FIELDS.get(tag, (tag, six.text_type))
        if rtype == list:
            d[tag] = [xml_to_dict(grandchild) for grandchild in child]
        elif rtype == dict:
            d[tag] = xml_to_dict(child)
        elif child.text is not None:
            d[tag] = rtype(child.text.strip())
    return d


class MassSpecApi(BaseChemSpider):

    def get_databases(self):
        """Get the list of datasources in ChemSpider."""
        response = self.request('MassSpecApi', 'GetDatabases')
        return [el.text for el in response]

    def get_extended_compound_info(self, csid):
        """Get extended record details for a CSID. Security token is required.

        :param string|int csid: ChemSpider ID.
        """
        response = self.request('MassSpecApi', 'GetExtendedCompoundInfo', csid=csid)
        return xml_to_dict(response)

    def get_extended_compound_info_list(self, csids):
        """Get extended record details for a list of CSIDs. Security token is required.

        :param list[string|int] csids: ChemSpider IDs.
        """
        response = self.request('MassSpecApi', 'GetExtendedCompoundInfoArray', csids=csids)
        return [xml_to_dict(result) for result in response]

    def get_extended_mol_compound_info_list(self, csids, mol_type=MOL2D, include_reference_counts=False,
                                            include_external_references=False):
        """Get extended record details (including MOL) for a list of CSIDs.

        A maximum of 250 CSIDs can be fetched per request. Security token is required.

        :param list[string|int] csids: ChemSpider IDs.
        :param string mol_type: :data:`~chemspipy.api.MOL2D`, :data:`~chemspipy.api.MOL3D` or
                                :data:`~chemspipy.api.BOTH`.
        :param bool include_reference_counts: Whether to include reference counts.
        :param bool include_external_references: Whether to include external references.
        """
        response = self.request('MassSpecApi', 'GetExtendedMolCompoundInfoArray', csids=csids,
                                eMolType=DIMENSIONS.get(mol_type, mol_type),
                                includeReferenceCounts=include_reference_counts,
                                includeExternalReferences=include_external_references)
        return [xml_to_dict(result) for result in response]

    def get_record_mol(self, csid, calc3d=False):
        """Get ChemSpider record in MOL format. Security token is required.

        :param string|int csid: ChemSpider ID.
        :param bool calc3d: Whether 3D coordinates should be calculated before returning record data.
        """
        response = self.request('MassSpecApi', 'GetRecordMol', csid=csid, calc3d=calc3d)
        return response.text

    def simple_search_by_formula(self, formula):
        """Search ChemSpider by molecular formula.

        :param string formula: Molecular formula
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        warnings.warn("Use search_by_formula instead of simple_search_by_formula.", DeprecationWarning)
        response = self.request('MassSpecApi', 'SearchByFormula2', formula=formula)
        return [Compound(self, el.text) for el in response]

    def simple_search_by_mass(self, mass, mass_range):
        """Search ChemSpider by mass +/- range.

        :param float mass: The mass to search for.
        :param float mass_range: The +/- mass range to allow.
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        warnings.warn("Use search_by_mass instead of simple_search_by_mass.", DeprecationWarning)
        response = self.request('MassSpecApi', 'SearchByMass2', mass=mass, range=mass_range)
        return [Compound(self, el.text) for el in response]

    # def get_compressed_records_sdf(self, rid):
    #     """Get an SDF containing all the results from a search operation.
    #
    #     A maximum of 10000 records can be fetched per request. Subscriber role security token is required.
    #
    #     Warning: This doesn't work reliably.
    #
    #     :param string rid: A transaction ID, returned by an asynchronous search method.
    #     :returns: SDF containing the requested records.
    #     :rtype: string
    #     """
    #     response = self.request('MassSpecApi', 'GetCompressedRecordsSdf', rid=rid, eComp='eGzip')
    #     if response.text:
    #         return zlib.decompress(b64decode(response.text.encode('utf-8')), 16+zlib.MAX_WBITS)
    #
    # def get_records_sdf(self, rid):
    #     """Get an SDF containing all the results from a search operation.
    #
    #     A maximum of 10000 records can be fetched per request. Subscriber role security token is required.
    #
    #     Warning: This doesn't work reliably.
    #
    #     :param string rid: A transaction ID, returned by an asynchronous search method.
    #     :returns: SDF containing the requested records.
    #     :rtype: string
    #     """
    #     response = self.request('MassSpecApi', 'GetRecordsSdf', rid=rid)
    #     if response.text:
    #         return response.text.encode('utf-8')


class SearchApi(BaseChemSpider):

    def async_simple_search(self, query):
        """Search ChemSpider with arbitrary query, returning results in order of the best match found.

        This method returns a transaction ID which can be used with other methods to get search status and results.

        Security token is required.

        :param string query: Search query - a name, SMILES, InChI, InChIKey, CSID, etc.
        :returns: Transaction ID.
        :rtype: string
        """
        response = self.request('Search', 'AsyncSimpleSearch', query=query)
        return response.text

    def async_simple_search_ordered(self, query, order=CSID, direction=ASCENDING):
        """Search ChemSpider with arbitrary query, returning results with a custom order.

        This method returns a transaction ID which can be used with other methods to get search status and results.

        Security token is required.

        :param string query: Search query - a name, SMILES, InChI, InChIKey, CSID, etc.
        :param string order: :data:`~chemspipy.api.CSID`, :data:`~chemspipy.api.MASS_DEFECT`,
                             :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
                             :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
                             :data:`~chemspipy.api.RSC_COUNT`.
        :param string direction: :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :returns: Transaction ID.
        :rtype: string
        """
        response = self.request('Search', 'AsyncSimpleSearchOrdered', query=query, orderBy=ORDERS[order],
                                orderDirection=DIRECTIONS[direction])
        return response.text

    def get_async_search_status(self, rid):
        """Check the status of an asynchronous search operation.

        Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :returns: Unknown, Created, Scheduled, Processing, Suspended, PartialResultReady, ResultReady, Failed,
                  TooManyRecords
        :rtype: string
        """
        response = self.request('Search', 'GetAsyncSearchStatus', rid=rid)
        return response.text

    def get_async_search_status_and_count(self, rid):
        """Check the status of an asynchronous search operation. If ready, a count and message are also returned.

        Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :rtype: dict
        """
        response = self.request('Search', 'GetAsyncSearchStatusAndCount', rid=rid)
        return xml_to_dict(response)

    def get_async_search_result(self, rid):
        """Get the results from a asynchronous search operation. Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        response = self.request('Search', 'GetAsyncSearchResult', rid=rid)
        return [Compound(self, el.text) for el in response]

    def get_async_search_result_part(self, rid, start=0, count=-1):
        """Get a slice of the results from a asynchronous search operation. Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :param int start: The number of results to skip.
        :param int count: The number of results to return. -1 returns all through to end.
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        response = self.request('Search', 'GetAsyncSearchResultPart', rid=rid, start=start, count=count)
        return [Compound(self, el.text) for el in response]

    def get_compound_info(self, csid):
        """Get SMILES, StdInChI and StdInChIKey for a given CSID. Security token is required.

        :param string|int csid: ChemSpider ID.
        :rtype: dict
        """
        response = self.request('Search', 'GetCompoundInfo', csid=csid)
        return xml_to_dict(response)

    def get_compound_thumbnail(self, csid):
        """Get PNG image as binary data.

        :param string|int csid: ChemSpider ID.
        :rtype: bytes
        """
        response = self.request('Search', 'GetCompoundThumbnail', id=csid)
        return b64decode(response.text.encode('utf-8'))

    def simple_search(self, query):
        """Search ChemSpider with arbitrary query.

        A maximum of 100 results are returned. Security token is required.

        :param string query: Search query - a name, SMILES, InChI, InChIKey, CSID, etc.
        :returns: List of :class:`Compounds <chemspipy.Compound>`.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        response = self.request('Search', 'SimpleSearch', query=query)
        return [Compound(self, el.text) for el in response]


class SpectraApi(BaseChemSpider):

    def get_all_spectra_info(self):
        """Get full list of all spectra in ChemSpider. Subscriber role security token is required.

        rtype: list[dict]
        """
        response = self.request('Spectra', 'GetAllSpectraInfo')
        return [xml_to_dict(result) for result in response]

    def get_spectrum_info(self, spectrum_id):
        """Get information for a specific spectrum ID. Subscriber role security token is required.

        :param string|int spectrum_id: spectrum ID.
        :returns: Spectrum info.
        :rtype: dict
        """
        response = self.request('Spectra', 'GetSpectrumInfo', spc_id=spectrum_id)
        return xml_to_dict(response)

    def get_compound_spectra_info(self, csid):
        """Get information about all the spectra for a ChemSpider ID. Subscriber role security token is required.

        :param string|int csid: ChemSpider ID.
        :returns: List of spectrum info.
        :rtype: list[dict]
        """
        response = self.request('Spectra', 'GetCompoundSpectraInfo', csid=csid)
        return [xml_to_dict(result) for result in response]

    def get_spectra_info_list(self, csids):
        """Get information about all the spectra for a list of ChemSpider IDs.

        :param list[string|int] csids: ChemSpider IDs.
        :returns: List of spectrum info.
        :rtype: list[dict]
        """
        response = self.request('Spectra', 'GetSpectraInfoArray', csids=csids)
        return [xml_to_dict(result) for result in response]


class InchiApi(BaseChemSpider):

    def get_original_mol(self, csid):
        """Get original submitted MOL file. Security token is required.

        :param string|int csid: ChemSpider ID.
        """
        response = self.request('InChI', 'CSIDToMol', csid=csid)
        return response.text

    def get_csid_from_inchikey(self, inchikey):
        """Get ChemSpiderID relating to an InChIKey.

        :param string inchi_key: InChIKey.
        """
        response = self.request('InChI', 'InChIKeyToCSID', inchi_key=inchikey)
        return response.text

    def get_mol_from_inchikey(self, inchikey):
        """Get mol relating to an InChIKey.

        :param string inchi_key: InChIKey.
        """
        response = self.request('InChI', 'InChIKeyToMol', inchi_key=inchikey)
        return response.text

    def get_inchi_from_inchikey(self, inchikey):
        """Get inchi by look-up of an InChIKey.

        :param string inchi_key: InChIKey.
        """
        response = self.request('InChI', 'InChIKeyToInChI', inchi_key=inchikey)
        return response.text

    # DONE
    # InChIKeyToCSID - inchi_key - csid
    # InChIKeyToMol - inchi_key - Mol
    # InChIKeyToInChI - inchi_key - InChI

    # TODO


    # InChIToCSID - inchi - csid
    # InChIToInChIKey - inchi - inchikey
    # InChIToMol - inchi - mol
    # InChIToSMILES - inchi - smiles
    # IsValidInChIKey - inchi_key - bool
    # MolToInChI - mol - inchi
    # MolToInChIKey - mol - inchi
    # ResolveInChIKey - inchi_key, out_format (MOL/SDF/SMILES/InChI) - list of strings
    # SMILESToInChI - smiles - inchi


class CustomApi(BaseChemSpider):

    def get_compound(self, csid):
        """Return a Compound object for a given ChemSpider ID. Security token is required.

        :param string|int csid: ChemSpider ID.
        :returns: The Compound with the specified ChemSpider ID.
        :rtype: :class:`~chemspipy.Compound`
        """
        return Compound(self, csid)

    def get_compounds(self, csids):
        """Return a list of Compound objects, given a list ChemSpider IDs. Security token is required.

        :param list[string|int] csids: List of ChemSpider IDs.
        :returns: List of Compounds with the specified ChemSpider IDs.
        :rtype: list[:class:`~chemspipy.Compound`]
        """
        return [Compound(self, csid) for csid in csids]

    def get_spectrum(self, spectrum_id):
        """Return a :class:`~chemspipy.Spectrum` object for a given spectrum ID. Subscriber role security token is required.

        :param string|int spectrum_id: Spectrum ID.
        :returns: The Spectrum with the specified spectrum ID.
        :rtype: :class:`~chemspipy.Spectrum`
        """
        return Spectrum(self, spectrum_id)

    def get_spectra(self, spectrum_ids):
        """Return a :class:`~chemspipy.Spectrum` object for a given spectrum ID. Subscriber role security token is required.

        :param list[string|int] spectrum_ids: List of spectrum IDs.
        :returns: List of spectra with the specified spectrum IDs.
        :rtype: list[:class:`~chemspipy.Spectrum`]
        """
        return [Spectrum(self, spectrum_id) for spectrum_id in spectrum_ids]

    def get_compound_spectra(self, csid):
        """Return :class:`~chemspipy.Spectrum` objects for all the spectra associated with a ChemSpider ID.

        :param csid: string|int csid: ChemSpider ID.
        :returns: List of spectra for the specified ChemSpider ID.
        :rtype: list[:class:`~chemspipy.Spectrum`]
        """
        return [Spectrum.from_info_dict(self, info) for info in self.get_spectra_info_list([csid])]

    def get_all_spectra(self):
        """Return a full list of :class:`~chemspipy.Spectrum` objects for all spectra in ChemSpider.

        Subscriber role security token is required.

        :returns: Full list of spectra in ChemSpider.
        :rtype: list[:class:`~chemspipy.Spectrum`]
        """
        return [Spectrum.from_info_dict(self, info) for info in self.get_all_spectra_info()]

    def search(self, query, order=None, direction=ASCENDING, raise_errors=False):
        """Search ChemSpider for the specified query and return the results. Security token is required.

        :param string|int query: Search query.
        :param string order: (Optional) :data:`~chemspipy.api.CSID`, :data:`~chemspipy.api.MASS_DEFECT`,
                             :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
                             :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
                             :data:`~chemspipy.api.RSC_COUNT`.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :param bool raise_errors: If True, raise exceptions. If False, store on Results ``exception`` property.
        :returns: Search Results list.
        :rtype: Results
        """
        if order and direction:
            return Results(self, self.async_simple_search_ordered, (query, order, direction), raise_errors=raise_errors)
        else:
            return Results(self, self.async_simple_search, (query,), raise_errors=raise_errors)

    # TODO: Wrappers for subscriber role asynchronous searches


class ChemSpider(CustomApi, MassSpecApi, SearchApi, SpectraApi, InchiApi):
    """Provides access to the ChemSpider API.

    Usage::

        >>> from chemspipy import ChemSpider
        >>> cs = ChemSpider('<YOUR-SECURITY-TOKEN>')

    """

    def __repr__(self):
        return 'ChemSpider()'
