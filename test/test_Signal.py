#!/usr/bin/env python
# Time-stamp: <2019-12-18 17:02:57 taoliu>

"""Module Description: Test functions for Signal.pyx

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest
import pytest

from math import log10
import numpy as np
import scipy
from MACS2.Signal import maxima, minima, savitzky_golay, savitzky_golay_order2_deriv1

# ------------------------------------
# Main function
# ------------------------------------

class Test_maxima(unittest.TestCase):

    def setUp(self):
        # the following signal array has two local maxima and one local minimum
        self.signal = np.array( [0] * 40 + [1] * 20 + [2] * 10 + [3] * 5 + [4] * 3 + [5] * 2 + [4] * 3 + [3] * 5 + [2] * 10 + [1] * 10 + [2] * 10 + [3] * 5 + [4] * 3 + [5] * 2 + [4] * 3 + [3] * 5 + [2] * 10 + [1] * 20 + [0] * 40, dtype="float32"  )
        self.windowsize = 41
        self.summit = [ 78, 126 ]
        self.trough = 102

    def test_maxima(self):
        expect = self.summit
        result = maxima( self.signal, self.windowsize )
        self.assertEqual( result[0], expect[0], msg=f"Not equal: result: {result}, expected: {expect}" )
        self.assertEqual( result[1], expect[1], msg=f"Not equal: result: {result}, expected: {expect}" )        

    def test_minima(self):
        expect = self.trough
        result = minima( self.signal, self.windowsize )[0]
        self.assertEqual( result, expect, msg=f"Not equal: result: {result}, expected: {expect}" )        

    def assertEqual_nparray1d ( self, a, b, places = 7 ):
        l = len(b)
        for i in range( l ):
            self.assertAlmostEqual( a[i], b[i], places = places, msg=f"Not equal at {i} {a[i]} {b[i]}" )
            
