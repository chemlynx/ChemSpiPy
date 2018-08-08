.. ChemSpiPy documentation master file, created by sphinx-quickstart on Mon Sep 15 10:51:47 2014.

ChemSpiPy
=========

.. sectionauthor:: Matt Swain <m.swain@me.com>

**ChemSpiPy** provides a way to interact with ChemSpider in Python. It allows chemical searches, chemical file
downloads, depiction and retrieval of chemical properties. Here's a quick peek::

    >>> from chemspipy import ChemSpider
    >>> cs = ChemSpider('<YOUR-API-KEY>')
    >>> c1 = cs.get_compound(236)  # Specify compound by ChemSpider ID
    >>> c2 = cs.search('benzene')  # Search using name, SMILES, InChI, InChIKey, etc.


Features
--------

- Search compounds by synonym, SMILES, InChI, InChIKey, formula and mass.
- Get identifiers and calculated properties for any compound record in ChemSpider.
- Download compound records as a MOL file with 2D or 3D coordinates.
- Get a 2D compound depiction as a PNG image.
- Complete interface to every endpoint of the ChemSpider Web APIs.
- Supports Python versions 2.7 – 3.4.

User guide
----------

A step-by-step guide to getting started with ChemSpiPy.

.. toctree::
   :maxdepth: 2

   guide/intro
   guide/install
   guide/gettingstarted
   guide/compound
   guide/searching
   guide/misc
   guide/advanced
   guide/contributing

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method.

.. toctree::
   :maxdepth: 2

   api
