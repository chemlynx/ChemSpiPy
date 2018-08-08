.. _api:

API documentation
=================

.. sectionauthor:: Matt Swain <m.swain@me.com>

.. module:: chemspipy

This part of the documentation is automatically generated from the ChemSpiPy source code and comments.

.. automodule:: chemspipy.api
.. autoclass:: chemspipy.ChemSpider()

   .. automethod:: get_compound(csid)
   .. automethod:: get_compounds(csids)
   .. automethod:: search(query, order=None, direction=ASCENDING, raise_errors=False)
   .. automethod:: simple_search(query)
   .. automethod:: get_record_mol(csid, calc3d=False)
   .. automethod:: get_compound_thumbnail(csid)
   .. automethod:: get_databases()
   .. automethod:: get_compound_info(csid)
   .. automethod:: get_extended_compound_info(csid)
   .. automethod:: get_extended_compound_info_list(csids)
   .. automethod:: get_extended_mol_compound_info_list(csids, mol_type=MOL2D, include_reference_counts=False, include_external_references=False)
   .. automethod:: request(api, endpoint, **params)
   .. automethod:: construct_api_url(api, endpoint, **params)
   .. automethod:: async_simple_search(query)
   .. automethod:: async_simple_search_ordered(query, order=CSID, direction=ASCENDING)
   .. automethod:: get_async_search_status(rid)
   .. automethod:: get_async_search_status_and_count(rid)
   .. automethod:: get_async_search_result(rid)
   .. automethod:: get_async_search_result_part(rid, start=0, count=-1)
   .. automethod:: simple_search_by_formula(formula)
   .. automethod:: simple_search_by_mass(mass, mass_range)

.. autodata:: chemspipy.api.MOL2D
.. autodata:: chemspipy.api.MOL3D
.. autodata:: chemspipy.api.BOTH
.. autodata:: chemspipy.api.ASCENDING
.. autodata:: chemspipy.api.DESCENDING
.. autodata:: chemspipy.api.CSID
.. autodata:: chemspipy.api.MASS_DEFECT
.. autodata:: chemspipy.api.MOLECULAR_WEIGHT
.. autodata:: chemspipy.api.REFERENCE_COUNT
.. autodata:: chemspipy.api.DATASOURCE_COUNT
.. autodata:: chemspipy.api.PUBMED_COUNT
.. autodata:: chemspipy.api.RSC_COUNT

.. automodule:: chemspipy.objects
.. autoclass:: chemspipy.Compound()
   :members:

.. automodule:: chemspipy.search
.. autoclass:: chemspipy.Results()
   :members:

.. automodule:: chemspipy.errors
.. autoexception:: chemspipy.errors.ChemSpiPyError()
.. autoexception:: chemspipy.errors.ChemSpiPyParseError()
.. autoexception:: chemspipy.errors.ChemSpiPyAuthError()
.. autoexception:: chemspipy.errors.ChemSpiPyNotFoundError()
.. autoexception:: chemspipy.errors.ChemSpiPyTimeoutError()
.. autoexception:: chemspipy.errors.ChemSpiPyServerError()
