.. _gettingstarted:

Getting started
===============

This page gives a introduction on how to get started with ChemSpiPy.

Before we start
---------------

- Make sure you have :ref:`installed ChemSpiPy <install>`.
- :ref:`Obtain an API key <apikey>` from the ChemSpider web site.

First steps
-----------

Start by importing ChemSpider::

    >>> from chemspipy import ChemSpider

Then connect to ChemSpider by creating a ``ChemSpider`` instance using your API key::

    >>> cs = ChemSpider('<YOUR-API-KEY>')

All your interaction with the ChemSpider database should now happen through this ChemSpider object, ``cs``.

Retrieve a Compound
-------------------

Retrieving information about a specific Compound in the ChemSpider database is simple.

Let's get the Compound with `ChemSpider ID 2157`_::

    >>> c = cs.get_compound(2157)

Now we have a :class:`~chemspipy.Compound` object called ``c``. We can get various identifiers and calculated
properties from this object::

    >>> print(c.molecular_formula)
    C_{9}H_{8}O_{4}
    >>> print(c.molecular_weight)
    180.15742
    >>> print(c.smiles)
    CC(=O)OC1=CC=CC=C1C(=O)O
    >>> print(c.common_name)
    Aspirin

Search for a name
-----------------

What if you don't know the ChemSpider ID of the Compound you want? Instead use the ``search`` method::

    >>> for result in cs.search('Glucose'):
    ...    print(result)
    Compound(5589)
    Compound(58238)
    Compound(71358)
    Compound(96749)
    Compound(9312824)
    Compound(9484839)

The ``search`` method accepts any identifer that ChemSpider can interpret, including names, registry numbers, SMILES
and InChI.

That's a quick taster of the basic ChemSpiPy functionality. Read on for more some more advanced usage examples.

.. _`ChemSpider ID 2157`: http://www.chemspider.com/Chemical-Structure.2157.html
