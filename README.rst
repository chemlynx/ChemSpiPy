ChemSpiPy
=========

.. image:: http://img.shields.io/pypi/v/ChemSpiPy.svg?style=flat
    :target: https://pypi.python.org/pypi/ChemSpiPy

.. image:: http://img.shields.io/pypi/l/ChemSpiPy.svg?style=flat
    :target: https://github.com/mcs07/ChemSpiPy/blob/master/LICENSE

.. image:: http://img.shields.io/travis/mcs07/ChemSpiPy/master.svg?style=flat
    :target: https://travis-ci.org/mcs07/ChemSpiPy

.. image:: http://img.shields.io/coveralls/mcs07/ChemSpiPy/master.svg?style=flat
    :target: https://coveralls.io/r/mcs07/ChemSpiPy?branch=master

ChemSpiPy provides a way to interact with ChemSpider in Python. It allows chemical searches, chemical file downloads,
depiction and retrieval of chemical properties::

    >>> from chemspipy import ChemSpider
    >>> cs = ChemSpider('<YOUR-API-KEY>')
    >>> c1 = cs.get_compound(236)  # Specify compound by ChemSpider ID
    >>> c2 = cs.search('benzene')  # Search using name, SMILES, InChI, InChIKey, etc.

Installation
------------

Install ChemSpiPy using pip::

    pip install chemspipy

Alternatively, try one of the other `installation options`_.

Documentation
-------------

Full documentation is available at http://chemspipy.readthedocs.io.

The `general documentation for the ChemSpider API`_ is also a useful resource.

Contribute
----------

-  Feature ideas and bug reports are welcome on the `Issue Tracker`_.
-  Fork the `source code`_ on GitHub, make changes and file a pull request.

License
-------

ChemSpiPy is licensed under the `MIT license`_.

This project was originally forked from `ChemSpiPy by Cameron Neylon`_, which has been released into the public domain.

.. _`installation options`: http://chemspipy.readthedocs.io/en/latest/guide/install.html
.. _`source code`: https://github.com/mcs07/ChemSpiPy
.. _`Issue Tracker`: https://github.com/mcs07/ChemSpiPy/issues
.. _`MIT license`: https://github.com/mcs07/ChemSpiPy/blob/master/LICENSE
.. _`ChemSpiPy by Cameron Neylon`: https://github.com/cameronneylon/ChemSpiPy
.. _`general documentation for the ChemSpider API`: http://www.chemspider.com/AboutServices.aspx
