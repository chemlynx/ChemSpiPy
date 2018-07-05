# -*- coding: utf-8 -*-
"""
chemspipy.objects
~~~~~~~~~~~~~~~~~

Objects returned by ChemSpiPy API methods.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

from .utils import memoized_property, timestamp


class Compound(object):
    """ A class for retrieving and caching details about a specific ChemSpider record.

    The purpose of this class is to provide access to various parts of the ChemSpider API that return information about
    a compound given its ChemSpider ID. Information is loaded lazily when requested, and cached for future access.
    """

    def __init__(self, cs, csid):
        """

        :param ChemSpider cs: ``ChemSpider`` session.
        :param int|string csid: ChemSpider ID.
        """
        self._cs = cs
        self._csid = int(csid)
        # TODO: Allow optional initialize  with a record-type response from the API (kwarg or class method from_dict?).

    def __eq__(self, other):
        return isinstance(other, Compound) and self.csid == other.csid

    def __repr__(self):
        return 'Compound(%r)' % self.csid

    def _repr_png_(self):
        """For IPython notebook, display 2D image."""
        return self.image

    @property
    def csid(self):
        """ChemSpider ID."""
        return self._csid

    # TODO: csid setter that clears cached properties?

    @property
    def image_url(self):
        """Return the URL of a PNG image of the 2D chemical structure."""
        return 'http://www.chemspider.com/ImagesHandler.ashx?id=%s' % self.csid

    @memoized_property
    def _compound_info(self):
        """Request compound info and cache the result."""
        return self._cs.get_compound_info(self.csid)

    @memoized_property
    def _extended_compound_info(self):
        """Request extended compound info and cache the result."""
        return self._cs.get_extended_compound_info(self.csid)

    @property
    def molecular_formula(self):
        """Return the molecular formula for this Compound.

        :rtype: string
        """
        return self._extended_compound_info['molecular_formula']

    @property
    def smiles(self):
        """Return the SMILES for this Compound.

        :rtype: string
        """
        return self._compound_info['smiles']

    @property
    def stdinchi(self):
        """Return the Standard InChI for this Compound.

        :rtype: string
        """
        return self._compound_info['inchi']

    @property
    def stdinchikey(self):
        """Return the Standard InChIKey for this Compound.

        :rtype: string
        """
        return self._compound_info['inchikey']

    @property
    def inchi(self):
        """Return the InChI for this Compound.

        :rtype: string
        """
        return self._extended_compound_info['inchi']

    @property
    def inchikey(self):
        """Return the InChIKey for this Compound.

        :rtype: string
        """
        return self._extended_compound_info['inchikey']

    @property
    def average_mass(self):
        """Return the average mass of this Compound.

        :rtype: float
        """
        return self._extended_compound_info['average_mass']

    @property
    def molecular_weight(self):
        """Return the molecular weight of this Compound.

        :rtype: float
        """
        return self._extended_compound_info['molecular_weight']

    @property
    def monoisotopic_mass(self):
        """Return the monoisotopic mass of this Compound.

        :rtype: float
        """
        return self._extended_compound_info['monoisotopic_mass']

    @property
    def nominal_mass(self):
        """Return the nominal mass of this Compound.

        :rtype: float
        """
        return self._extended_compound_info['nominal_mass']

    @property
    def alogp(self):
        """Return the calculated AlogP for this Compound.

        :rtype: float
        """
        return self._extended_compound_info['alogp']

    @property
    def xlogp(self):
        """Return the calculated XlogP for this Compound.

        :rtype: float
        """
        return self._extended_compound_info['xlogp']

    @property
    def common_name(self):
        """Return the common name for this Compound.

        :rtype: string
        """
        return self._extended_compound_info['common_name']

    @memoized_property
    def mol_2d(self):
        """Return the MOL file for this Compound with 2D coordinates.

        :rtype: string
        """
        return self._cs.get_record_mol(self.csid, calc3d=False)

    @memoized_property
    def mol_3d(self):
        """Return the MOL file for this Compound with 3D coordinates.

        :rtype: string
        """
        return self._cs.get_record_mol(self.csid, calc3d=True)

    @memoized_property
    def mol_raw(self):
        """Return unprocessed MOL file for this Compound.

        :rtype: string
        """
        return self._cs.get_original_mol(self.csid)

    @memoized_property
    def image(self):
        """Return a 2D depiction of this Compound.

        :rtype: bytes
        """
        return self._cs.get_compound_thumbnail(self.csid)

'''
    @memoized_property
    def spectra(self):
        """Return all the available spectral data for this Compound.

        :rtype: list[:class:`~chemspipy.Spectrum`]
        """
        return [Spectrum.from_info_dict(self._cs, info) for info in self._cs.get_spectra_info_list([self.csid])]



class Spectrum(object):
    """ A class for retrieving and caching details about a Spectrum."""

    def __init__(self, cs, spectrum_id):
        """Initializing a Spectrum from a spectrum ID requires a subscriber role security token.

        :param ChemSpider cs: ``ChemSpider`` session.
        :param int|string spectrum_id: Spectrum ID.
        """
        self._cs = cs
        self._spectrum_id = int(spectrum_id)

    def __eq__(self, other):
        return isinstance(other, Spectrum) and self.spectrum_id == other.spectrum_id

    def __repr__(self):
        return 'Spectrum(%r)' % self.spectrum_id

    @classmethod
    def from_info_dict(cls, cs, info):
        """Initialize a Spectrum from an info dict that has already been retrieved."""
        s = cls(cs, info['spectrum_id'])
        s._info = info
        return s

    @property
    def _spectrum_info(self):
        """Full spectrum info.

        :rtype: dict
        """
        if not hasattr(self, '_info'):
            self._info = self._cs.get_spectrum_info(self._spectrum_id)
        return self._info

    @property
    def spectrum_id(self):
        """Spectrum ID.

        :rtype: int
        """
        return self._spectrum_id

    @property
    def csid(self):
        """ChemSpider ID of related compound.

        :rtype: int
        """
        return self._spectrum_info['csid']

    @property
    def spectrum_type(self):
        """Spectrum type.

        Possible values include HNMR, CNMR, IR, UV-Vis, NIR, EI, 2D1H1HCOSY, 2D1H13CD, APCI+, R, MALDI+, 2D1H13CLR,
        APPI-, CI+ve, ESI+, 2D1H1HOESY, FNMR, CI-ve, ESI-, PNMR.

        :rtype: string
        """
        return self._spectrum_info['spectrum_type']

    @property
    def file_name(self):
        """Spectrum file name.

        :rtype: string
        """
        return self._spectrum_info['file_name']

    @property
    def comments(self):
        """Spectrum comments. Can be None.

        :rtype: string
        """
        return self._spectrum_info.get('comments')

    @property
    def url(self):
        """Spectrum URL.

        :rtype: string
        """
        return 'https://www.chemspider.com/FilesHandler.ashx?type=blob&disp=1&id=%s' % self.spectrum_id

    @memoized_property
    def data(self):
        """Spectrum data file contents. Requires an additional request. Result is cached.

        :rtype: string
        """
        r = self._cs.http.get(self.url)
        return r.text

    @property
    def original_url(self):
        """Original spectrum URL. Can be None.

        :rtype: string
        """
        return self._spectrum_info.get('original_url')

    @property
    def submitted_date(self):
        """Spectrum submitted date.

        :rtype: :py:class:`datetime.datetime`
        """
        return timestamp(self._spectrum_info['submitted_date'])
'''
