.. _advanced:

Advanced
========

Keep your API key secret
------------------------

Be careful not to include your API key when sharing code. A simple way to ensure this doesn't happen by accident is to
store your API key as an environment variable that can be specified in your `.bash_profile` or `.zshrc` file::

    export CHEMSPIDER_API_KEY=<YOUR-API-KEY>

This can then be retrieved in your scripts using ``os.environ``::

    >>> api_key = os.environ['CHEMSPIDER_API_KEY']
    >>> cs = ChemSpider(api_key)

Specify a User Agent
--------------------

As well as using your API key, it is possible to identify your program to the ChemSpider servers using a User
Agent string.

You can specify a custom User Agent through ChemSpiPy through the optional ``user_agent`` parameter to the ChemSpider
class::

    >>> from chemspipy import ChemSpider
    >>> cs = ChemSpider('<YOUR-API-KEY>', user_agent='My program 1.3, ChemSpiPy 1.0.5, Python 3.6')

Logging
-------

ChemSpiPy can generate logging statements if required. Just set the desired logging level::

    import logging
    logging.basicConfig(level=logging.DEBUG)

The logger is named 'chemspipy'. There is more information on logging in the `Python logging documentation`_.

.. _`Python logging documentation`: http://docs.python.org/2/howto/logging.html
