# -*- coding: utf-8 -*-
"""
chemspipy.api
~~~~~~~~~~~~~

Core API for interacting with ChemSpider web services.

:copyright: Copyright 2018 by Matt Swain, David Sharpe.
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

try:
    import json
except ImportError:
    print('json module not loaded')

import requests
import six

from . import __version__
from .errors import ChemSpiPyError, ChemSpiPyParseError, ChemSpiPyAuthError, ChemSpiPyServerError
from .errors import ChemSpiPyNotFoundError
from .objects import Compound #, Spectrum
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

    def __init__(self, apikey=None, user_agent=None, api_base_url='https://api.rsc.org', api_path='compounds', version='v1'):
        """

        :param string security_token: (Optional) Your ChemSpider security token.
        :param string apikey: (Required) Your developer portal apikey.
        :param string user_agent: (Optional) Identify your application to ChemSpider servers.
        :param string api_url: (Optional) Alternative API server.
        """
        log.debug('Initializing ChemSpider')
        self.api_url = '%s/%s/%s' % (api_base_url, api_path, version)
        self.http = requests.session()
        self.http.headers['User-Agent'] = user_agent if user_agent else 'ChemSpiPy/%s Python/%s ' % (__version__, sys.version.split()[0])
        self.apikey = apikey
        self.http.headers['Content-Type'] = 'application/json'

    def request(self, api, endpoint, http_verb, record_id=None, **params):
        """Construct API request and return the json response.

        :param string api: The specific ChemSpider API to call (MassSpec, Search, Spectra, InChI).
        :param string new_api: The specific ChemSpider API to call (compounds).
        :param string new_api: The version of the ChemSpider API to call (v1).
        :param string endpoint: ChemSpider API endpoint.
        :param params: (Optional) Parameters for the ChemSpider endpoint as keyword arguments.
        :rtype: xml tree, json
        """
        # "https://api.rsc.org/compounds/v1/lookups/datasources"
        if record_id:
            url = '%s/%s/%s/%s' % (self.api_url, api, record_id, endpoint)
            if 'fields' in params.keys():
                print(params['fields'])
                parameters = ','.join(params['fields'])
                url = '%s?fields=%s' % (url, parameters)
        else:

            url = '%s/%s/%s' % (self.api_url, api, endpoint)

        print(url)
        log.debug('Request: %s %s', url, params)
        params['apikey'] = self.apikey

        if http_verb == 'POST':
            url = '%s/%s/%s' % (self.api_url, api, endpoint)
            self.http.headers.update({'apikey': self.apikey})
            try:
                print(params['data'])
                response = self.http.post(url, json=params['data'], headers=self.http.headers)
            except requests.RequestException as e:
                raise ChemSpiPyError(six.text_type(e))

        elif http_verb == 'GET':
            #url = construct_api_url('MassSpecAPI', 'GetExtendedCompoundInfo', csid=0)
            self.http.headers.update({'apikey': self.apikey})
            try:
                response = requests.get(url, headers=self.http.headers)
            except requests.RequestException as e:
                raise ChemSpiPyError(six.text_type(e))

        if response.status_code == 200:
            pass
        elif response.status_code == 400:
            raise ChemSpiPyAuthError('Bad Request. Check the request you sent and try again.')
        elif response.status_code == 401:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Unauthorized. Check you have supplied the correct API key and that you have sent it as an HTTP Header called 'apikey'.")
        elif response.status_code == 404:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Not Found. The requested endpoint URL is not recognized. Change your request and try again.")
        elif response.status_code == 405:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Method Not Allowed. The verb is incorrect for the endpoint. Change your request and try again.")
        elif response.status_code == 413:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Payload Too Large. The request you sent was too big to handle. Change your request and try again.")
        elif response.status_code == 429:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Too Many Requests. Send fewer requests, or use rate-limiting to slow them down, then try again.")
        elif response.status_code == 500:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Internal Server Error. Wait and try again.")
        elif response.status_code == 503:
            # Generally when supplying a security token with incorrect format
            raise ChemSpiPyAuthError("Service Unavailable. Wait and try again.")
        else:
            raise ChemSpiPyServerError(response.text)

        try:
            response_dict = json.loads(response.text)
        except json.JSONDecodeError as e:
            raise ChemSpiPyParseError('Unable to parse json response: %s' % e)

        return response_dict



class MassSpecApi(BaseChemSpider):

    def get_databases(self):
        """Get the list of datasources in ChemSpider."""
        response = self.request( 'lookups', 'datasources', 'GET')


        response_list = response['dataSources']
        print(response_list)
        return response_list

    def get_extended_compound_info(self, csid, **params):
        """Get extended record details for a CSID. Security token is required.

        :param string|int csid: ChemSpider ID.
        :param list[string] fields: types of data to return, mol formula etc.
        """
        response = self.request('records', 'details', 'GET', record_id=csid, fields=params['fields'])

        print(response)
        return response

    '''
    def get_extended_compound_info_list(self, csids):
        """Get extended record details for a list of CSIDs. Security token is required.

        :param list[string|int] csids: ChemSpider IDs.
        """
        response = self.request('records', 'batch', 'POST' csids=csids)
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
    '''

    def get_record_mol(self, csid, calc3d=False):
        """Get ChemSpider record in MOL format. Security token is required.

        :param string|int csid: ChemSpider ID.
        :param bool calc3d: Whether 3D coordinates should be calculated before returning record data.
        """
        if calc3d==False:
            response = self.request('records', 'details', 'GET', record_id=csid, fields=['mol2D'])
        else:
            response = self.request('records', 'details', 'GET', record_id=csid, fields=['mol3D'])

        keys = response.keys()
        for key in keys:
            if key != 'id':
                break

        output = response.pop(key)
        return output

    def simple_search_by_formula(self, formula, datasources=None):
        """Search ChemSpider by molecular formula.

        :param string formula: Molecular formula
        :param list datasources: datasources
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """

        payload = {
                "formula": formula
                }
        if datasources:
            data['dataSources'] = datasources

        response = self.request('filter', 'formula', 'POST', data=payload)
        print(response)
        return [Compound(self, el.text) for el in response]

    def simple_search_by_mass(self, mass, mass_range, datasources=None):
        """Search ChemSpider by mass +/- range.

        :param float mass: The mass to search for.
        :param float mass_range: The +/- mass range to allow.
        :param list datasources: datasources
        :returns: A list of Compounds.
        :rtype: list[:class:`~chemspipy.Compound`]
        """

        payload = {
                "mass": mass,
                "range": mass_range
                }
        if datasources:
            data['dataSources'] = datasources

        response = self.request('filter', 'mass', 'POST', data=payload)

        print(response)
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

        payload = {
                "name": query
                }

        response = self.request('filter', 'name', 'POST', data=payload)
        return response

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
        payload = {
                "name": query,
                'orderBy': order,
                'orderDirection': direction
                }


        response = self.request('filter', 'name', 'POST', data=payload)
        return response

    def get_async_search_status(self, rid):
        """Check the status of an asynchronous search operation.

        Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :returns: Unknown, Created, Scheduled, Processing, Suspended, PartialResultReady, ResultReady, Failed,
                  TooManyRecords
        :rtype: string
        """
        response = self.request('filter', 'status', 'GET', record_id=rid)
        return response

    def get_async_search_status_and_count(self, rid):
        """Check the status of an asynchronous search operation. If ready, a count and message are also returned.

        Security token is required.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :rtype: dict
        """
        response = self.request('filter', 'status', 'GET', record_id=rid)
        return response

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

    def get_csid_from_inchi(self, inchi):
        """Get ChemSpider ID by look-up of an InChI.

        :param string inchi: InChI.
        """
        response = self.request('InChI', 'InChIToCSID', inchi=inchi)
        return response.text

    def get_inchikey_from_inchi(self, inchi):
        """Get InChIKey by look-up of an InChI.

        :param string inchi: InChI.
        """
        response = self.request('InChI', 'InChIToInChIKey', inchi=inchi)
        return response.text

    def get_mol_from_inchi(self, inchi):
        """Get mol by look-up of an InChI.

        :param string inchi: InChI.
        """
        response = self.request('InChI', 'InChIToMol', inchi=inchi)
        return response.text

    # DONE
    # InChIKeyToCSID - inchi_key - csid
    # InChIKeyToMol - inchi_key - Mol
    # InChIKeyToInChI - inchi_key - InChI
    # InChIToCSID - inchi - csid
    # InChIToInChIKey - inchi - inchikey
    # InChIToMol - inchi - mol

    # TODO

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


class ChemSpider(CustomApi, MassSpecApi, SearchApi, InchiApi):
    """Provides access to the ChemSpider API.

    Usage::

        >>> from chemspipy import ChemSpider
        >>> cs = ChemSpider('<YOUR-SECURITY-TOKEN>')

    """

    def __repr__(self):
        return 'ChemSpider()'
