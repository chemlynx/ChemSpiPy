.. _compound:

Compound
========

Many ChemSpiPy search methods return :class:`~chemspipy.Compound` objects, which provide more functionality that a
simple list of ChemSpider IDs. The primary benefit is allowing easy access to further compound properties after
performing a search.

Creating a Compound
-------------------

The easiest way to create a :class:`~chemspipy.Compound` for a given ChemSpider ID is to use the ``get_compound``
method::

    >>> compound = cs.get_compound(2157)

Alternatively, a :class:`~chemspipy.Compound` can be instantiated directly::

    >>> compound = Compound(cs, 2157)

Either way, no requests are made to the ChemSpider servers until specific :class:`~chemspipy.Compound` properties are
requested::

    >>> print(compound.molecular_formula)
    C_{9}H_{8}O_{4}
    >>> print(compound.molecular_weight)
    180.15742
    >>> print(compound.smiles)
    CC(=O)OC1=CC=CC=C1C(=O)O
    >>> print(compound.common_name)
    Aspirin

Properties are cached locally after the first time they are retrieved, speeding up subsequent access and reducing the
number of unnecessary requests to the ChemSpider servers.

Searching for Compounds
-----------------------

See the :ref:`searching documentation <searching>` for full details.

Compound properties
-------------------

- ``csid``: ChemSpider ID.
- ``image_url``: URL of a PNG image of the 2D chemical structure.
- ``molecular_formula``: Molecular formula.
- ``smiles``: SMILES string.
- ``stdinchi``: Standard InChI string.
- ``stdinchikey``: Standard InChIKey.
- ``inchi``: InChI string.
- ``inchikey``: InChIKey.
- ``average_mass``: Average mass.
- ``molecular_weight``: Molecular weight.
- ``monoisotopic_mass``: Monoisotopic mass.
- ``nominal_mass``: Nominal mass.
- ``alogp``: AlogP.
- ``xlogp``: XlogP.
- ``common_name``: Common Name.
- ``mol_2d``: MOL file containing 2D coordinates.
- ``mol_3d``: MOL file containing 3D coordinates.
- ``mol_raw``: Unprocessed MOL file.
- ``image``: 2D depiction as binary data in PNG format.

Implementation details
----------------------

Each :class:`~chemspipy.Compound` object is a simple wrapper around a ChemSpider ID. Behind the scenes, the property
methods use the ``get_compound_info``, ``get_extended_compound_info``, ``get_record_mol`` and
``get_compound_thumbnail`` API methods to retrieve the relevant information. It is possible to use these API methods
directly if required::

    >>> info = cs.get_extended_compound_info(2157)
    {u'smiles': u'CC(=O)Oc1ccccc1C(=O)O', u'common_name': u'Aspirin', u'nominal_mass': 180.0, u'molecular_formula': u'C_{9}H_{8}O_{4}', u'inchikey': u'BSYNRYMUTXBXSQ-UHFFFAOYAW', u'molecular_weight': 180.1574, u'inchi': u'InChI=1/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)', u'average_mass': 180.1574, u'csid': 2157, u'alogp': 0.0, u'xlogp': 0.0, u'monoisotopic_mass': 180.042252}

Results are returned as a python dictionary that is derived directly from the ChemSpider API XML response.
