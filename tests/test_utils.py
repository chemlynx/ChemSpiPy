# -*- coding: utf-8 -*-
"""
test_utils
~~~~~~~~~~

Test miscellaneous utility functions.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import datetime
import logging

from chemspipy.utils import timestamp, duration


logging.basicConfig(level=logging.WARN, format='%(levelname)s:%(name)s:(%(threadName)-10s):%(message)s')
logging.getLogger('chemspipy').setLevel(logging.DEBUG)


def test_timestamp_microseconds():
    """Test timestamp parser function on timestamps strings with microseconds."""
    assert timestamp('2007-08-08T20:18:36.593') == datetime.datetime(2007, 8, 8, 20, 18, 36, 593000)
    assert timestamp('2035-12-31T23:59:59.999') == datetime.datetime(2035, 12, 31, 23, 59, 59, 999000)


def test_timestamp_seconds():
    """Test timestamp parser function on timestamps strings with no microseconds."""
    assert timestamp('2007-08-08T20:18:36') == datetime.datetime(2007, 8, 8, 20, 18, 36)
    assert timestamp('2010-09-01T04:33:59') == datetime.datetime(2010, 9, 1, 4, 33, 59)


def test_duration_microseconds():
    """Test duration parser function on duration strings with microseconds."""
    assert duration('0:00:00.001') == datetime.timedelta(0, 0, 1000)
    assert duration('0:00:00.12') == datetime.timedelta(0, 0, 120000)
    assert duration('0:00:00.120') == datetime.timedelta(0, 0, 120000)
    assert duration('0:00:00.052') == datetime.timedelta(0, 0, 52000)
    assert duration('0:00:03.523') == datetime.timedelta(0, 3, 523000)


def test_duration_seconds():
    """Test duration parser function on duration strings with no microseconds."""
    assert duration('0:00:00') == datetime.timedelta(0)
    assert duration('0:00:03') == datetime.timedelta(0, 3)
