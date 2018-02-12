from unittest import TestCase

from atmopy import compute_vars

class testcomputevars(TestCase):
    def test_compute_PR(self):
        pr,atts = compute_vars.compute_PR('./atmopy/tests/wrfout_test.nc')
        self.assertTrue(pr.mean()>=0.)
